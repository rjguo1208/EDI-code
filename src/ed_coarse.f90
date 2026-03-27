MODULE ed_coarse
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ed_coarse_calc, load_supercell_pot, read_filukk_edi
  PUBLIC :: ed_fine_interp, ed_fine_interp_offdiag, read_edmatw_file
  PUBLIC :: read_edmatw_2d_file
  PUBLIC :: ed_interp_from_file
  PUBLIC :: ed_coarse_full_q
  PUBLIC :: ed_fine_interp_2d
  PUBLIC :: ed_direct_from_files, generate_nscf_input

CONTAINS

  SUBROUTINE ed_coarse_calc(nbndsub, nrr_k, irvec_k, ndegen_k, wslen_k, &
                             edmatw, prefix_in)
    !-----------------------------------------------------------------------
    ! Part A: Compute electron-defect matrix elements on the coarse k-grid.
    !
    ! Potentials V_d, V_p, V_colin must be pre-loaded in edic_mod
    ! before calling this routine. The primitive cell state (dffts, igk_k,
    ! ngk, xk, nks) must be current (from read_file).
    !
    ! Workflow:
    !   1. Allocate workspace (evc1, evc2, psic1, psic2 in edic_mod)
    !   2. Loop k: read_collected_wfc, calcmdefect_ml_rs + calcmdefect_mnl_ks
    !   3. Bloch -> Wannier via edbloch2wane
    !   4. Write edmatw and decay.M
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast, mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_global, ONLY : inter_pool_comm
    USE cell_base, ONLY : at, bg, alat
    USE klist, ONLY : xk, igk_k, ngk, nkstot, nks
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module, ONLY : npol
    USE fft_base, ONLY : dffts
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE wann_common, ONLY : u_mat, u_mat_opt, n_wannier, &
                             num_bands, excluded_band
    USE ep_constants, ONLY : czero, cone, zero, twopi, ci, bohr2ang
    USE edbloch2wan, ONLY : edbloch2wane
    USE edic_mod, ONLY : V_d, V_p, V_colin, &
                          evc1, evc2, psic1, psic2
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr_k
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k), ndegen_k(nrr_k)
    REAL(dp), INTENT(IN) :: wslen_k(nrr_k)
    COMPLEX(dp), INTENT(OUT) :: edmatw(nbndsub, nbndsub, nrr_k)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, ibnd, jbnd, ib_i, ib_j
    INTEGER :: nbnd_kept
    COMPLEX(dp) :: mlocal, mnonlocal0, mnonlocal1
    COMPLEX(dp), ALLOCATABLE :: edmatkq(:,:,:), cu(:,:,:)
    REAL(dp), ALLOCATABLE :: xk_loc(:,:)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part A: Computing e-d matrix on coarse grid'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I6)') 'Supercell FFT grid: ', &
            V_d%nr1, V_d%nr2, V_d%nr3
       WRITE(stdout, '(5X,A,I6)') 'Supercell nat (defect): ', V_d%nat
       WRITE(stdout, '(5X,A,I6)') 'Primitive nkstot = ', nkstot
       WRITE(stdout, '(5X,A,I6)') 'Primitive nks    = ', nks
    ENDIF

    nbnd_kept = 0
    DO ibnd = 1, nbnd
       IF (.NOT. excluded_band(ibnd)) nbnd_kept = nbnd_kept + 1
    ENDDO

    ! Allocate edic_mod workspace (used by calcmdefect_ml_rs/mnl_ks)
    ALLOCATE(evc1(npwx * npol, nbnd))
    ALLOCATE(evc2(npwx * npol, nbnd))
    ALLOCATE(psic1(dffts%nnr))
    ALLOCATE(psic2(dffts%nnr))

    ALLOCATE(edmatkq(nbnd_kept, nbnd_kept, nks))
    ALLOCATE(cu(nbnd_kept, nbndsub, nks))
    ALLOCATE(xk_loc(3, nks))
    xk_loc(:, 1:nks) = xk(:, 1:nks)

    ! Build combined rotation matrix: cu = U_opt * U_dis
    DO ik = 1, nks
       IF (num_bands > n_wannier) THEN
          cu(:,:,ik) = MATMUL(u_mat_opt(:,:,ik), u_mat(:,:,ik))
       ELSE
          cu(:,:,ik) = u_mat(:,:,ik)
       ENDIF
    ENDDO

    edmatkq = czero

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,I6,A)') &
            'Computing M for ', nks, ' local k-points...'
       WRITE(stdout, '(5X,A,I4,A,I4,A)') &
            '  (', nbnd_kept, ' x ', nbnd_kept, ' band pairs per k-point)'
       FLUSH(stdout)
    ENDIF

    CALL compute_edmatkq_fast(nbnd, nbnd_kept, nks, edmatkq)

    !===========================================================
    ! Bloch -> Wannier transform
    !===========================================================
    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'Rotating M to Wannier gauge and FT to real space...'
       FLUSH(stdout)
    ENDIF

    CALL edbloch2wane(nbnd_kept, nbndsub, nks, nkstot, xk_loc, cu, cu, &
                       edmatkq, nrr_k, irvec_k, wslen_k, edmatw)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'Writing edmatw to file...'
       CALL write_edmatw(nbndsub, nrr_k, irvec_k, ndegen_k, edmatw, prefix_in)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

    DEALLOCATE(edmatkq, cu, xk_loc)
    DEALLOCATE(evc1, evc2, psic1, psic2)

  END SUBROUTINE ed_coarse_calc

  SUBROUTINE compute_edmatkq_fast(nbnd, nbnd_kept, nks, edmatkq)
    !-----------------------------------------------------------------------
    ! Optimized matrix element computation for coarse grid (k = k').
    !
    ! Key optimizations:
    !   1. Fold supercell potential onto primitive cell grid ONCE
    !      (17.3M → 80K points, ~200× fewer per band pair)
    !   2. FFT each band to real space ONCE per k-point
    !   3. Local M = dot product over primitive grid (not supercell loop)
    !   4. Nonlocal: compute beta projectors ONCE per k-point
    !-----------------------------------------------------------------------
    USE kinds, ONLY : dp
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE wvfct, ONLY : npwx
    USE noncollin_module, ONLY : npol, noncolin, lspinorb
    USE klist, ONLY : xk, igk_k, ngk
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE wann_common, ONLY : excluded_band
    USE uspp_param, ONLY : nh
    USE uspp, ONLY : dvan, dvan_so
    USE becmod, ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
    USE edic_mod, ONLY : V_file, V_d, V_p, V_colin, evc1, evc2
    USE ep_constants, ONLY : czero
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbnd, nbnd_kept, nks
    COMPLEX(dp), INTENT(OUT) :: edmatkq(nbnd_kept, nbnd_kept, nks)

    INTEGER :: ik, ibnd, ib_i, ib_j, ig, irx, iry, irz, ir, ipol
    INTEGER :: ir1mod, ir2mod, ir3mod, inr, irnmod
    INTEGER :: na, nt, ih, jh, ikb, jkb, ijkb0, nkb_d, nkb_p
    COMPLEX(dp) :: mlocal, mnl_d, mnl_p
    COMPLEX(dp), ALLOCATABLE :: psir(:,:,:)  ! (nnr, npol, nbnd_kept)
    COMPLEX(dp), ALLOCATABLE :: psic_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:)
    REAL(dp), ALLOCATABLE :: V_folded(:)
    TYPE(bec_type) :: becp_d, becp_p
    INTEGER, ALLOCATABLE :: band_map(:)

    ! Build band mapping
    ALLOCATE(band_map(nbnd_kept))
    ib_i = 0
    DO ibnd = 1, nbnd
       IF (excluded_band(ibnd)) CYCLE
       ib_i = ib_i + 1
       band_map(ib_i) = ibnd
    ENDDO

    ! ============================================================
    ! FOLD supercell potential onto primitive cell grid (done ONCE)
    ! V_folded(r_prim) = sum_{R_super -> r_prim} V_colin(R_super)
    ! ============================================================
    ALLOCATE(V_folded(dffts%nnr))
    V_folded(:) = 0.0_dp
    inr = 0
    DO irz = 0, V_d%nr3 - 1
       ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
       DO iry = 0, V_d%nr2 - 1
          ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
          DO irx = 0, V_d%nr1 - 1
             ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
             inr = inr + 1
             irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
             V_folded(irnmod) = V_folded(irnmod) + V_colin(inr)
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,I8,A,I8)') &
            'V_folded: supercell grid=', V_d%nr1*V_d%nr2*V_d%nr3, &
            ' -> primitive grid=', dffts%nnr
       FLUSH(stdout)
    ENDIF

    ! Count nonlocal projectors for defect and pristine supercells
    nkb_d = 0
    DO nt = 1, V_d%ntyp
       DO na = 1, V_d%nat
          IF (V_d%ityp(na) == nt) nkb_d = nkb_d + nh(nt)
       ENDDO
    ENDDO
    nkb_p = 0
    DO nt = 1, V_p%ntyp
       DO na = 1, V_p%nat
          IF (V_p%ityp(na) == nt) nkb_p = nkb_p + nh(nt)
       ENDDO
    ENDDO

    ALLOCATE(psir(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psic_tmp(dffts%nnr))
    ALLOCATE(vkb_d(npwx, nkb_d))
    ALLOCATE(vkb_p(npwx, nkb_p))

    DO ik = 1, nks
       ! Read wavefunctions
       CALL read_collected_wfc(restart_dir(), ik, evc1)
       evc2(:,:) = evc1(:,:)

       ! FFT all kept bands to real space ONCE
       ! For SOC (npol=2): FFT both spinor components separately
       DO ib_i = 1, nbnd_kept
          ibnd = band_map(ib_i)
          DO ipol = 1, npol
             psic_tmp(:) = (0.0_dp, 0.0_dp)
             DO ig = 1, ngk(ik)
                psic_tmp(dffts%nl(igk_k(ig, ik))) = evc1(ig + (ipol-1)*npwx, ibnd)
             ENDDO
             CALL invfft('Wave', psic_tmp, dffts)
             psir(:, ipol, ib_i) = psic_tmp(:)
          ENDDO
       ENDDO

       ! Compute beta projectors ONCE per k-point
       CALL get_betavkb(ngk(ik), igk_k(1,ik), xk(1,ik), &
                         vkb_d, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
       CALL allocate_bec_type(nkb_d, nbnd, becp_d)
       CALL calbec(ngk(ik), vkb_d, evc1, becp_d)

       CALL get_betavkb(ngk(ik), igk_k(1,ik), xk(1,ik), &
                         vkb_p, V_p%nat, V_p%ityp, V_p%tau, nkb_p)
       CALL allocate_bec_type(nkb_p, nbnd, becp_p)
       CALL calbec(ngk(ik), vkb_p, evc1, becp_p)

       ! Compute matrix elements using folded potential + pre-computed bec
       DO ib_i = 1, nbnd_kept
          DO ib_j = 1, nbnd_kept
             ! LOCAL: sum over spinor components
             ! M_local = Σ_σ Σ_r ψ*_{σ}(r) ΔV(r) ψ'_{σ}(r)
             mlocal = czero
             DO ipol = 1, npol
                DO ir = 1, dffts%nnr
                   mlocal = mlocal + CONJG(psir(ir, ipol, ib_i)) * &
                        psir(ir, ipol, ib_j) * V_folded(ir)
                ENDDO
             ENDDO
             mlocal = mlocal / DBLE(dffts%nnr)

             ! NONLOCAL defect: use pre-computed becp_d
             mnl_d = czero
             ijkb0 = 0
             DO nt = 1, V_d%ntyp
                DO na = 1, V_d%nat
                   IF (V_d%ityp(na) == nt) THEN
                      IF (lspinorb) THEN
                         ! SOC: 4 spin channels using dvan_so and becp%nc
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            mnl_d = mnl_d + &
                                 CONJG(becp_d%nc(ikb,1,band_map(ib_i))) * becp_d%nc(ikb,1,band_map(ib_j)) * dvan_so(ih,ih,1,nt) + &
                                 CONJG(becp_d%nc(ikb,1,band_map(ib_i))) * becp_d%nc(ikb,2,band_map(ib_j)) * dvan_so(ih,ih,2,nt) + &
                                 CONJG(becp_d%nc(ikb,2,band_map(ib_i))) * becp_d%nc(ikb,1,band_map(ib_j)) * dvan_so(ih,ih,3,nt) + &
                                 CONJG(becp_d%nc(ikb,2,band_map(ib_i))) * becp_d%nc(ikb,2,band_map(ib_j)) * dvan_so(ih,ih,4,nt)
                            DO jh = ih + 1, nh(nt)
                               jkb = ijkb0 + jh
                               mnl_d = mnl_d + &
                                    CONJG(becp_d%nc(ikb,1,band_map(ib_i))) * becp_d%nc(jkb,1,band_map(ib_j)) * dvan_so(ih,jh,1,nt) + &
                                    CONJG(becp_d%nc(jkb,1,band_map(ib_i))) * becp_d%nc(ikb,1,band_map(ib_j)) * dvan_so(jh,ih,1,nt) + &
                                    CONJG(becp_d%nc(ikb,1,band_map(ib_i))) * becp_d%nc(jkb,2,band_map(ib_j)) * dvan_so(ih,jh,2,nt) + &
                                    CONJG(becp_d%nc(jkb,1,band_map(ib_i))) * becp_d%nc(ikb,2,band_map(ib_j)) * dvan_so(jh,ih,2,nt) + &
                                    CONJG(becp_d%nc(ikb,2,band_map(ib_i))) * becp_d%nc(jkb,1,band_map(ib_j)) * dvan_so(ih,jh,3,nt) + &
                                    CONJG(becp_d%nc(jkb,2,band_map(ib_i))) * becp_d%nc(ikb,1,band_map(ib_j)) * dvan_so(jh,ih,3,nt) + &
                                    CONJG(becp_d%nc(ikb,2,band_map(ib_i))) * becp_d%nc(jkb,2,band_map(ib_j)) * dvan_so(ih,jh,4,nt) + &
                                    CONJG(becp_d%nc(jkb,2,band_map(ib_i))) * becp_d%nc(ikb,2,band_map(ib_j)) * dvan_so(jh,ih,4,nt)
                            ENDDO
                         ENDDO
                      ELSE
                         ! Scalar-relativistic: becp%k and dvan
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            mnl_d = mnl_d + CONJG(becp_d%k(ikb, band_map(ib_i))) * &
                                 becp_d%k(ikb, band_map(ib_j)) * dvan(ih, ih, nt)
                            DO jh = ih + 1, nh(nt)
                               jkb = ijkb0 + jh
                               mnl_d = mnl_d + &
                                    (CONJG(becp_d%k(ikb, band_map(ib_i))) * becp_d%k(jkb, band_map(ib_j)) + &
                                     CONJG(becp_d%k(jkb, band_map(ib_i))) * becp_d%k(ikb, band_map(ib_j))) * &
                                    dvan(ih, jh, nt)
                            ENDDO
                         ENDDO
                      ENDIF
                      ijkb0 = ijkb0 + nh(nt)
                   ENDIF
                ENDDO
             ENDDO

             ! NONLOCAL pristine: use pre-computed becp_p
             mnl_p = czero
             ijkb0 = 0
             DO nt = 1, V_p%ntyp
                DO na = 1, V_p%nat
                   IF (V_p%ityp(na) == nt) THEN
                      IF (lspinorb) THEN
                         ! SOC: 4 spin channels using dvan_so and becp%nc
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            mnl_p = mnl_p + &
                                 CONJG(becp_p%nc(ikb,1,band_map(ib_i))) * becp_p%nc(ikb,1,band_map(ib_j)) * dvan_so(ih,ih,1,nt) + &
                                 CONJG(becp_p%nc(ikb,1,band_map(ib_i))) * becp_p%nc(ikb,2,band_map(ib_j)) * dvan_so(ih,ih,2,nt) + &
                                 CONJG(becp_p%nc(ikb,2,band_map(ib_i))) * becp_p%nc(ikb,1,band_map(ib_j)) * dvan_so(ih,ih,3,nt) + &
                                 CONJG(becp_p%nc(ikb,2,band_map(ib_i))) * becp_p%nc(ikb,2,band_map(ib_j)) * dvan_so(ih,ih,4,nt)
                            DO jh = ih + 1, nh(nt)
                               jkb = ijkb0 + jh
                               mnl_p = mnl_p + &
                                    CONJG(becp_p%nc(ikb,1,band_map(ib_i))) * becp_p%nc(jkb,1,band_map(ib_j)) * dvan_so(ih,jh,1,nt) + &
                                    CONJG(becp_p%nc(jkb,1,band_map(ib_i))) * becp_p%nc(ikb,1,band_map(ib_j)) * dvan_so(jh,ih,1,nt) + &
                                    CONJG(becp_p%nc(ikb,1,band_map(ib_i))) * becp_p%nc(jkb,2,band_map(ib_j)) * dvan_so(ih,jh,2,nt) + &
                                    CONJG(becp_p%nc(jkb,1,band_map(ib_i))) * becp_p%nc(ikb,2,band_map(ib_j)) * dvan_so(jh,ih,2,nt) + &
                                    CONJG(becp_p%nc(ikb,2,band_map(ib_i))) * becp_p%nc(jkb,1,band_map(ib_j)) * dvan_so(ih,jh,3,nt) + &
                                    CONJG(becp_p%nc(jkb,2,band_map(ib_i))) * becp_p%nc(ikb,1,band_map(ib_j)) * dvan_so(jh,ih,3,nt) + &
                                    CONJG(becp_p%nc(ikb,2,band_map(ib_i))) * becp_p%nc(jkb,2,band_map(ib_j)) * dvan_so(ih,jh,4,nt) + &
                                    CONJG(becp_p%nc(jkb,2,band_map(ib_i))) * becp_p%nc(ikb,2,band_map(ib_j)) * dvan_so(jh,ih,4,nt)
                            ENDDO
                         ENDDO
                      ELSE
                         ! Scalar-relativistic: becp%k and dvan
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            mnl_p = mnl_p + CONJG(becp_p%k(ikb, band_map(ib_i))) * &
                                 becp_p%k(ikb, band_map(ib_j)) * dvan(ih, ih, nt)
                            DO jh = ih + 1, nh(nt)
                               jkb = ijkb0 + jh
                               mnl_p = mnl_p + &
                                    (CONJG(becp_p%k(ikb, band_map(ib_i))) * becp_p%k(jkb, band_map(ib_j)) + &
                                     CONJG(becp_p%k(jkb, band_map(ib_i))) * becp_p%k(ikb, band_map(ib_j))) * &
                                    dvan(ih, jh, nt)
                            ENDDO
                         ENDDO
                      ENDIF
                      ijkb0 = ijkb0 + nh(nt)
                   ENDIF
                ENDDO
             ENDDO

             edmatkq(ib_i, ib_j, ik) = mlocal + mnl_d - mnl_p
          ENDDO
       ENDDO

       CALL deallocate_bec_type(becp_d)
       CALL deallocate_bec_type(becp_p)

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,I6,A,I6,A,I6,A)') &
               '  k-point ', ik, ' / ', nks, ' done (', &
               nbnd_kept*nbnd_kept, ' matrix elements)'
          WRITE(stdout, '(5X,A,2ES14.6,A,2ES14.6)') &
               '    M(1,1)=', REAL(edmatkq(1,1,ik)), AIMAG(edmatkq(1,1,ik)), &
               '  M(1,2)=', REAL(edmatkq(1,MIN(2,nbnd_kept),ik)), &
               AIMAG(edmatkq(1,MIN(2,nbnd_kept),ik))
          FLUSH(stdout)
       ENDIF
    ENDDO

    DEALLOCATE(psir, psic_tmp, V_folded, vkb_d, vkb_p, band_map)

  END SUBROUTINE compute_edmatkq_fast

  SUBROUTINE ed_fine_interp(nbndsub, nrr, irvec, ndegen, chw, edmatw, &
                             nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B: Wannier → Bloch interpolation on fine k-grid.
    !
    ! For each k_fine on the fine grid:
    !   1. cfac = exp(ik_fine·R) for all R
    !   2. Diagonalize H(k_fine) → eigenvalues ε(k_fine), eigenvectors U(k_fine)
    !   3. Interpolate M: M_W(k_fine) = Σ_R exp(ik·R)/ndegen(R) · M(R)
    !   4. Rotate to Bloch: M_B = U† · M_W · U
    !   5. Output |M_B(n,m,k)|² for transport
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE constants, ONLY : rytoev
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)    ! H(R)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)  ! M(R)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, iunit
    REAL(dp) :: xk_cryst(3)
    REAL(dp) :: eig(nbndsub)
    COMPLEX(dp) :: evec(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)   ! M in Wannier basis at k
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub)   ! M in Bloch basis at k
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname

    nktotf = nk1f * nk2f * nk3f

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B: Wannier interpolation on fine grid'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Open output file for interpolated matrix elements
    iunit = 86
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmatf.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Interpolated e-d matrix elements |M(n,m,k)|^2'
       WRITE(iunit, '(A)') '# ik  kx  ky  kz  ibnd  jbnd  |M|^2(Ry^2)  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_cryst(1) = DBLE(i1) / DBLE(nk1f)
             xk_cryst(2) = DBLE(i2) / DBLE(nk2f)
             xk_cryst(3) = DBLE(i3) / DBLE(nk3f)

             ! Phase factors exp(ik·R)
             CALL get_cfac(nrr, irvec, xk_cryst, cfac)

             ! Diagonalize H(k) → eigenvalues + eigenvectors U(k)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig, evec, chw, cfac)

             ! Interpolate M(R) → M_W(k) in Wannier basis
             CALL edmatwan2bloch(nbndsub, nrr, ndegen, edmatw, cfac, edmatf_w)

             ! Rotate to Bloch basis: M_B = U† · M_W · U
             ! tmp = M_W · U
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec, nbndsub, czero, tmp, nbndsub)
             ! M_B = U† · tmp
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             ! Write to file
             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_cryst, &
                           ibnd, jbnd, &
                           ABS(edmatf_b(ibnd, jbnd))**2, &
                           REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,I8,A)') 'Interpolated M at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp

  SUBROUTINE ed_fine_interp_offdiag(nbndsub, nrr, irvec, ndegen, chw, edmatw, &
                                     nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B extension: Off-diagonal k-point interpolation.
    !
    ! For a fixed initial k_i, compute M(n,m, k_i, k_f) for ALL k_f
    ! on the fine grid. Uses momentum transfer q = k_f - k_i:
    !
    !   M_W(i,j,q) = Σ_R exp(iq·R)/ndegen(R) · M(i,j,R)
    !   M_B(n,m,k_i,k_f) = U†(k_i) · M_W(q) · U(k_f)
    !
    ! Output: prefix_edmat_scatter.dat
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, ir, iunit
    REAL(dp) :: xk_i(3), xk_f(3), xq(3)
    REAL(dp) :: eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_i(nrr), cfac_f(nrr), cfac_q(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub)
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname
    REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
    COMPLEX(dp), PARAMETER :: ci = (0.0_dp, 1.0_dp)

    ! K point in crystal coordinates (MoS2): (1/3, 1/3, 0)
    xk_i = (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp/)

    nktotf = nk1f * nk2f * nk3f

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B: Off-diagonal scattering from K point'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3F10.5)') 'Initial k (crystal): ', xk_i
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Diagonalize H(k_i) once → eigenvalues + eigenvectors at initial k
    CALL get_cfac(nrr, irvec, xk_i, cfac_i)
    CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_i, evec_i, chw, cfac_i)

    ! Find CBM index at K (band with smallest eigenvalue above gap)
    ! For MoS2 with 11 bands, CBM at K is typically band 8
    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'Eigenvalues at K (Ry):'
       DO ibnd = 1, nbndsub
          WRITE(stdout, '(5X,A,I3,A,F12.6)') '  band ', ibnd, ': ', eig_i(ibnd)
       ENDDO
       FLUSH(stdout)
    ENDIF

    ! Open output file
    iunit = 85
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmat_scatter.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Off-diagonal scattering: M(n, m, k_i=K, k_f)'
       WRITE(iunit, '(A,3F10.5)') '# k_initial (crystal) = ', xk_i
       WRITE(iunit, '(A)') '# ik  kfx  kfy  kfz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_f(1) = DBLE(i1) / DBLE(nk1f)
             xk_f(2) = DBLE(i2) / DBLE(nk2f)
             xk_f(3) = DBLE(i3) / DBLE(nk3f)

             ! Momentum transfer q = k_f - k_i
             xq = xk_f - xk_i

             ! Phase factors for q: exp(iq·R)
             DO ir = 1, nrr
                cfac_q(ir) = EXP(ci * twopi * DOT_PRODUCT(xq, DBLE(irvec(:, ir))))
             ENDDO

             ! Interpolate M(R) at q: M_W(q) = Σ_R exp(iq·R)/ndegen · M(R)
             CALL edmatwan2bloch(nbndsub, nrr, ndegen, edmatw, cfac_q, edmatf_w)

             ! Diagonalize H(k_f) → eigenvectors U(k_f)
             CALL get_cfac(nrr, irvec, xk_f, cfac_f)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_f)

             ! Rotate to Bloch: M_B = U†(k_i) · M_W(q) · U(k_f)
             ! tmp = M_W · U(k_f)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                          cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             ! M_B = U†(k_i) · tmp
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                          cone, evec_i, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             ! Write all band pairs
             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, &
                           ABS(edmatf_b(ibnd, jbnd))**2, &
                           REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,I8,A)') 'Computed M(K -> k_f) at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp_offdiag


  SUBROUTINE ed_coarse_full_q(nbndsub, nrr_k, irvec_k, ndegen_k, &
                               edmatw_2d, prefix_in)
    !-----------------------------------------------------------------------
    ! Full double-FT: M(k_i, k_f) for ALL (k_i, k_f) pairs on the coarse grid.
    ! Computes M(R_e, R_p) via (paper Eq.5):
    !   M_W(k_i, k_f) = U†(k_i) · M_B(k_i, k_f) · U(k_f)
    !   M(R_e, R_p) = (1/Nk²) Σ_{ki,kf} exp(+iki·R_e) exp(-ikf·R_p) M_W
    !
    ! Pool-parallel: each pool handles its local k_i, broadcasts ψ(k_f).
    ! No clean_pw, no read_file_new — QE state preserved.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_global, ONLY : inter_pool_comm
    USE mp_pools, ONLY : npool, my_pool_id, intra_pool_comm
    USE cell_base, ONLY : at
    USE klist, ONLY : xk, igk_k, ngk, nkstot, nks
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module, ONLY : npol
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE wann_common, ONLY : u_mat, n_wannier, num_bands, excluded_band
    USE edic_mod, ONLY : V_d, V_p, V_colin
    USE ep_constants, ONLY : czero, cone
    USE constants, ONLY : tpi
    USE parallelism, ONLY : fkbounds
    USE wigner_seitz_edi, ONLY : wigner_seitz_k
    USE parallel_include
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr_k
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k), ndegen_k(nrr_k)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw_2d(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, ib_i, ib_j, ir, ir_e, ir_p
    INTEGER :: irx, iry, irz, ir1mod, ir2mod, ir3mod, inr, irnmod
    INTEGER :: ibnd, nbnd_kept, lower_bnd, upper_bnd
    INTEGER :: ikf_p, ikf_l, ierr_mpi
    INTEGER :: my_rank, num_ranks, nkf, nrr_p
    REAL(dp) :: xkf_cryst(3), rdotk, rdotkf, arg, d1, d2, d3
    REAL(dp), ALLOCATABLE :: xk_cryst_loc(:,:), xk_cryst_all(:,:)
    COMPLEX(dp) :: mlocal, phase
    COMPLEX(dp), ALLOCATABLE :: V_folded_kf(:)
    COMPLEX(dp), ALLOCATABLE :: psir_k(:,:,:), psir_kq(:,:,:), psic_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: evc_tmp(:,:)
    COMPLEX(dp), ALLOCATABLE :: edmatkq(:,:), edms(:,:), eptmp(:,:)
    COMPLEX(dp), ALLOCATABLE :: cu_k(:,:), cu_kq(:,:)
    COMPLEX(dp), ALLOCATABLE :: cu_all(:,:,:)
    INTEGER, ALLOCATABLE :: band_map(:)

    ! q-grid = k-grid (every k-point is a valid momentum transfer)
    nkf = nkstot
    nrr_p = nrr_k  ! R_e and R_p on the same WS grid
    CALL MPI_Comm_rank(world_comm, my_rank, ierr_mpi)
    CALL MPI_Comm_size(world_comm, num_ranks, ierr_mpi)
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Full double-FT: M(k, k'') for all (k,k'') pairs'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6)') 'Nk_i = Nk_f = ', nkstot
       WRITE(stdout, '(5X,A,I4,A,I4)') 'nks=', nks, '  npool=', npool
       WRITE(stdout, '(5X,A,I6)') 'nrr (R_e = R_p)=', nrr_k
       FLUSH(stdout)
    ENDIF

    ! Build band mapping
    nbnd_kept = 0
    ALLOCATE(band_map(nbnd))
    DO ibnd = 1, nbnd
       IF (.NOT. excluded_band(ibnd)) THEN
          nbnd_kept = nbnd_kept + 1
          band_map(nbnd_kept) = ibnd
       ENDIF
    ENDDO

    ! Gather xk in crystal coords for ALL k-points (small array)
    ALLOCATE(xk_cryst_loc(3, nks), xk_cryst_all(3, nkstot))
    xk_cryst_loc(:, 1:nks) = xk(:, 1:nks)
    CALL cryst_to_cart(nks, xk_cryst_loc, at, -1)
    xk_cryst_all = 0.0_dp
    xk_cryst_all(:, lower_bnd:upper_bnd) = xk_cryst_loc(:, 1:nks)
    CALL mp_sum(xk_cryst_all, inter_pool_comm)
    CALL cryst_to_cart(nks, xk_cryst_loc, at, 1)  ! restore to Cartesian

    ! Gather cu (combined rotation) for ALL k-points (small array)
    ALLOCATE(cu_all(nbnd_kept, nbndsub, nkstot))
    cu_all = czero
    cu_all(:, :, lower_bnd:upper_bnd) = u_mat(:, :, 1:nks)
    CALL mp_sum(cu_all, inter_pool_comm)

    ! Cache pool-local wavefunctions in real space ONCE (avoid repeated disk I/O)
    ALLOCATE(psir_k(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psir_kq(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psic_tmp(dffts%nnr))
    ALLOCATE(evc_tmp(npwx * npol, nbnd))
    ALLOCATE(V_folded_kf(dffts%nnr))
    ALLOCATE(edmatkq(nbnd_kept, nbnd_kept))
    ALLOCATE(cu_k(nbnd_kept, nbndsub), cu_kq(nbnd_kept, nbndsub))
    ALLOCATE(edms(nbndsub, nbndsub), eptmp(nbnd_kept, nbndsub))
    ALLOCATE(edmatw_2d(nbndsub, nbndsub, nrr_k, nrr_p))
    edmatw_2d = czero

    ! Panel broadcast: cache local psi(k), compute V_folded on-the-fly per (ki,kf)
    BLOCK
       COMPLEX(dp), ALLOCATABLE :: psir_cache(:,:,:,:), psir_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: phase_ki(:), phase_kf(:)
       INTEGER :: ik_cache, ib_cache, ig_cache, ik_global, ipol_cache
       INTEGER :: ip, src_lower, src_nks, nkbase, nkrest, nks_max
       INTEGER :: ikf_local, ikf_global, npairs_done
       REAL(dp) :: qcryst(3), xki_cryst(3)

       ! Maximum k-points per pool (for receive buffer sizing)
       nkbase = nkstot / npool
       nkrest = MOD(nkstot, npool)
       nks_max = nkbase
       IF (nkrest > 0) nks_max = nkbase + 1

       ! Step 1: Cache pool-local wavefunctions in real space
       ! For SOC (npol=2): cache both spinor components
       ALLOCATE(psir_cache(dffts%nnr, npol, nbnd_kept, nks))
       DO ik_cache = 1, nks
          CALL read_collected_wfc(restart_dir(), ik_cache, evc_tmp)
          DO ib_cache = 1, nbnd_kept
             DO ipol_cache = 1, npol
                psic_tmp = (0.0_dp, 0.0_dp)
                DO ig_cache = 1, ngk(ik_cache)
                   psic_tmp(dffts%nl(igk_k(ig_cache, ik_cache))) = &
                        evc_tmp(ig_cache + (ipol_cache-1)*npwx, band_map(ib_cache))
                ENDDO
                CALL invfft('Wave', psic_tmp, dffts)
                psir_cache(:, ipol_cache, ib_cache, ik_cache) = psic_tmp
             ENDDO
          ENDDO
       ENDDO

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,I4,A)') &
               'Cached ', nks, ' pool-local k-points in real space'
          FLUSH(stdout)
       ENDIF

       ! Step 2: Panel broadcast — iterate over source pools in lockstep
       ! All pools call MPI_Bcast with the SAME root, preventing deadlock.
       ! V_folded is computed ON-THE-FLY for each (ki,kf) pair using the
       ! EXACT q = kf - ki (no BZ wrapping) to avoid exp(iG*r) phase errors.
       ALLOCATE(psir_recv(dffts%nnr, npol, nbnd_kept, nks_max))
       ALLOCATE(phase_ki(nrr_k), phase_kf(nrr_k))
       edmatw_2d = czero
       npairs_done = 0

       DO ip = 0, npool - 1
          ! Determine k-point range of source pool ip
          IF (ip < nkrest) THEN
             src_lower = ip * (nkbase + 1) + 1
             src_nks = nkbase + 1
          ELSE
             src_lower = nkrest * (nkbase + 1) + (ip - nkrest) * nkbase + 1
             src_nks = nkbase
          ENDIF

          ! Source pool broadcasts its cached wavefunctions to all pools
          IF (my_pool_id == ip) THEN
             psir_recv(:, :, :, 1:src_nks) = psir_cache(:, :, :, 1:nks)
          ENDIF
          CALL MPI_Bcast(psir_recv, dffts%nnr * npol * nbnd_kept * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)

          ! Process all (ki_local, kf_from_src) pairs
          DO ikf_local = 1, src_nks
             ikf_global = src_lower + ikf_local - 1

             DO ik = 1, nks
                ik_global = ik + lower_bnd - 1

                ! Compute V_folded on-the-fly with EXACT q = kf - ki (no wrapping!)
                ! Using wrapped q would give V_folded(q+G) = exp(iG*r) * V_folded(q),
                ! introducing an unphysical phase that corrupts the double-FT.
                qcryst = xk_cryst_all(:, ikf_global) - xk_cryst_all(:, ik_global)
                V_folded_kf = (0.0_dp, 0.0_dp)
                IF (ABS(qcryst(1)) < 1.0d-8 .AND. ABS(qcryst(2)) < 1.0d-8 &
                     .AND. ABS(qcryst(3)) < 1.0d-8) THEN
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + &
                                     ir2mod * dffts%nr1 + ir1mod + 1
                            V_folded_kf(irnmod) = V_folded_kf(irnmod) + V_colin(inr)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   ! Phase per grid step: q·r = 2π × (q1·irx/nr1 + q2·iry/nr2 + q3·irz/nr3)
                   ! q is in crystal coords, grid indices are in crystal-lattice units
                   ! DO NOT mix with at (Cartesian lattice vectors) — that gives wrong
                   ! phases for non-orthogonal lattices (hexagonal, etc.)
                   d1 = tpi * qcryst(1) / dffts%nr1
                   d2 = tpi * qcryst(2) / dffts%nr2
                   d3 = tpi * qcryst(3) / dffts%nr3
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + &
                                     ir2mod * dffts%nr1 + ir1mod + 1
                            arg = irx * d1 + iry * d2 + irz * d3
                            phase = CMPLX(COS(arg), SIN(arg), dp)
                            V_folded_kf(irnmod) = V_folded_kf(irnmod) + V_colin(inr) * phase
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

                ! Matrix element: M_B = sum_{r,σ} u*_{σ}(ki,r) u_{σ}(kf,r) V_q(r) / nnr
                psir_k(:,:,:) = psir_cache(:,:,:, ik)
                psir_kq(:,:,:) = psir_recv(:,:,:, ikf_local)

                DO ib_i = 1, nbnd_kept
                   DO ib_j = 1, nbnd_kept
                      mlocal = czero
                      DO ipol_cache = 1, npol
                         DO ir = 1, dffts%nnr
                            mlocal = mlocal + CONJG(psir_k(ir, ipol_cache, ib_i)) * &
                                 psir_kq(ir, ipol_cache, ib_j) * V_folded_kf(ir)
                         ENDDO
                      ENDDO
                      edmatkq(ib_i, ib_j) = mlocal / DBLE(dffts%nnr)
                   ENDDO
                ENDDO

                ! Wannier rotation: M_W = U_dag(ki) M_B U(kf)
                cu_k  = cu_all(:, :, ik_global)
                cu_kq = cu_all(:, :, ikf_global)
                CALL ZGEMM('N', 'N', nbnd_kept, nbndsub, nbnd_kept, &
                     cone, edmatkq, nbnd_kept, cu_kq, nbnd_kept, czero, eptmp, nbnd_kept)
                CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbnd_kept, &
                     cone, cu_k, nbnd_kept, eptmp, nbnd_kept, czero, edms, nbndsub)

                ! Double FT: M(Re,Rp) += exp(+iki*Re) exp(-ikf*Rp) / Nk^2 * M_W
                xki_cryst = xk_cryst_all(:, ik_global)
                xkf_cryst = xk_cryst_all(:, ikf_global)
                DO ir = 1, nrr_k
                   rdotk  = tpi * DOT_PRODUCT(xki_cryst, DBLE(irvec_k(:, ir)))
                   phase_ki(ir) = EXP(CMPLX(0.0_dp, +rdotk, dp))
                   rdotkf = tpi * DOT_PRODUCT(xkf_cryst, DBLE(irvec_k(:, ir)))
                   phase_kf(ir) = EXP(CMPLX(0.0_dp, -rdotkf, dp))
                ENDDO
                DO ir_p = 1, nrr_k
                   DO ir_e = 1, nrr_k
                      edmatw_2d(:,:,ir_e,ir_p) = edmatw_2d(:,:,ir_e,ir_p) + &
                           phase_ki(ir_e) * phase_kf(ir_p) / DBLE(nkstot * nkstot) * edms(:,:)
                   ENDDO
                ENDDO

                npairs_done = npairs_done + 1
             ENDDO  ! ik (local ki)
          ENDDO  ! ikf_local (kf from source pool)

          IF (ionode) THEN
             WRITE(stdout, '(5X,A,I4,A,I4,A,I6,A,I6)') &
                  'kf pool ', ip, ' / ', npool - 1, &
                  '  pairs: ', npairs_done, ' / ', nks * nkstot
             FLUSH(stdout)
          ENDIF
       ENDDO  ! ip (source pool)

       DEALLOCATE(phase_ki, phase_kf, psir_recv, psir_cache)

       ! Gather FT contributions from all pools
       CALL mp_sum(edmatw_2d, inter_pool_comm)
    END BLOCK



    DEALLOCATE(psir_k, psir_kq, psic_tmp, evc_tmp, V_folded_kf)
    DEALLOCATE(edmatkq, cu_k, cu_kq, edms, eptmp)
    DEALLOCATE(xk_cryst_loc, xk_cryst_all, cu_all, band_map)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'All ki-kf pairs completed. Writing output...'
       CALL write_edmatw_2d_v2(nbndsub, nrr_k, nrr_p, irvec_k, irvec_k, &
                                ndegen_k, ndegen_k, edmatw_2d, prefix_in)
       ! Write full M(R_e,R_p) binary for edwread
       CALL write_edmatw_2d_bin(nbndsub, nrr_k, irvec_k, ndegen_k, &
                                 edmatw_2d, prefix_in)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

  CONTAINS

    FUNCTION find_q_index(qvec, nk, xk_c) RESULT(iq_out)
      ! Find k-grid index closest to qvec (modulo reciprocal lattice)
      INTEGER, INTENT(IN) :: nk
      REAL(dp), INTENT(IN) :: qvec(3), xk_c(3, nk)
      INTEGER :: iq_out, jk
      REAL(dp) :: diff(3), dist, min_dist
      min_dist = 1.0d10
      iq_out = 1
      DO jk = 1, nk
         diff = qvec - xk_c(:, jk)
         diff = diff - NINT(diff)
         dist = SUM(diff**2)
         IF (dist < min_dist) THEN
            min_dist = dist
            iq_out = jk
         ENDIF
      ENDDO
    END FUNCTION

    FUNCTION find_kf_cryst(iki, dkf, nk, xk_c) RESULT(ikf_out)
      INTEGER, INTENT(IN) :: iki, nk
      REAL(dp), INTENT(IN) :: dkf(3), xk_c(3, nk)
      INTEGER :: ikf_out, jk
      REAL(dp) :: xkf_target(3), diff(3), dist, min_dist
      xkf_target = xk_c(:, iki) + dkf
      xkf_target = xkf_target - FLOOR(xkf_target)
      min_dist = 1.0d10
      ikf_out = 1
      DO jk = 1, nk
         diff = xkf_target - xk_c(:, jk)
         diff = diff - NINT(diff)
         dist = SUM(diff**2)
         IF (dist < min_dist) THEN
            min_dist = dist
            ikf_out = jk
         ENDIF
      ENDDO
    END FUNCTION

    FUNCTION pool_of_k(ik_global, nktot, np) RESULT(ipool)
      INTEGER, INTENT(IN) :: ik_global, nktot, np
      INTEGER :: ipool, nk_per_pool, remainder
      nk_per_pool = nktot / np
      remainder = MOD(nktot, np)
      IF (ik_global <= (nk_per_pool + 1) * remainder) THEN
         ipool = (ik_global - 1) / (nk_per_pool + 1)
      ELSE
         ipool = remainder + (ik_global - (nk_per_pool + 1) * remainder - 1) / nk_per_pool
      ENDIF
    END FUNCTION

    FUNCTION local_k_index(ik_global, nktot, np, ipool) RESULT(ik_local)
      INTEGER, INTENT(IN) :: ik_global, nktot, np, ipool
      INTEGER :: ik_local, nk_per_pool, remainder, offset
      nk_per_pool = nktot / np
      remainder = MOD(nktot, np)
      IF (ipool < remainder) THEN
         offset = ipool * (nk_per_pool + 1)
      ELSE
         offset = remainder * (nk_per_pool + 1) + (ipool - remainder) * nk_per_pool
      ENDIF
      ik_local = ik_global - offset
    END FUNCTION

    SUBROUTINE fetch_wfc_from_pool(src_pool, ik_src, evc_out, npwx_npol, nb, &
                                    pool_comm, my_pid, np)
      INTEGER, INTENT(IN) :: src_pool, ik_src, npwx_npol, nb, pool_comm, my_pid, np
      COMPLEX(dp), INTENT(OUT) :: evc_out(npwx_npol, nb)
      INTEGER :: tag, ierr, status(MPI_STATUS_SIZE)

      tag = ik_src + 10000
      IF (my_pid == src_pool) THEN
         ! I own this k-point — read and send to all who need it
         CALL read_collected_wfc(restart_dir(), ik_src, evc_out)
      ENDIF
      ! Broadcast from src_pool to all pools
      CALL MPI_Bcast(evc_out, npwx_npol * nb * 2, MPI_DOUBLE_PRECISION, &
                      src_pool, pool_comm, ierr)
    END SUBROUTINE

  END SUBROUTINE ed_coarse_full_q


  SUBROUTINE write_edmatw_2d(nbndsub, nrr, irvec, ndegen, edmatw_2d, prefix_in)
    USE io_global, ONLY : ionode, stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir_e, ir_p, ii, jj, iunit
    REAL(dp) :: max_abs
    CHARACTER(LEN=256) :: fname

    iunit = 81
    fname = TRIM(prefix_in) // '_edmatw_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# M(R_e, R_p) from full double-FT'
    WRITE(iunit, '(I6, I6)') nbndsub, nrr

    DO ir_p = 1, nrr
       DO ir_e = 1, nrr
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          WRITE(iunit, '(6I5, ES16.8)') irvec(:,ir_e), irvec(:,ir_p), max_abs
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'M(R_e, R_p) written to ', TRIM(fname)

    ! Also write decay comparison
    fname = TRIM(prefix_in) // '_decay_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# |R_e|  |R_p|  max|M(R_e,R_p)|  (for comparison with diagonal)'
    DO ir_p = 1, nrr
       DO ir_e = 1, nrr
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          IF (max_abs > 1.0d-10) THEN
             WRITE(iunit, '(2F14.6, ES16.8)') &
                  SQRT(DBLE(SUM(irvec(:,ir_e)**2))), &
                  SQRT(DBLE(SUM(irvec(:,ir_p)**2))), max_abs
          ENDIF
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'Decay data written to ', TRIM(fname)

  END SUBROUTINE write_edmatw_2d

  SUBROUTINE write_edmatw_2d_v2(nbndsub, nrr_e, nrr_p, irvec_e, irvec_p, &
                                 ndegen_e, ndegen_p, edmatw_2d, prefix_in)
    USE io_global, ONLY : ionode, stdout
    USE cell_base, ONLY : at, alat
    USE ep_constants, ONLY : bohr2ang
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr_e, nrr_p
    INTEGER, INTENT(IN) :: irvec_e(3, nrr_e), irvec_p(3, nrr_p)
    INTEGER, INTENT(IN) :: ndegen_e(nrr_e), ndegen_p(nrr_p)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr_e, nrr_p)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir_e, ir_p, iunit
    REAL(dp) :: max_abs, re_len, rp_len, rvec_e(3), rvec_p(3)
    CHARACTER(LEN=256) :: fname

    iunit = 81

    ! Write decay data (|R_e| in Angstrom, |R_p| in Angstrom, max|M|)
    fname = TRIM(prefix_in) // '_decay_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# |R_e|(Ang)  |R_p|(Ang)  max|M(R_e,R_p)|'
    DO ir_p = 1, nrr_p
       rvec_p = MATMUL(at, DBLE(irvec_p(:, ir_p))) * alat * bohr2ang
       rp_len = SQRT(DOT_PRODUCT(rvec_p, rvec_p))
       DO ir_e = 1, nrr_e
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          IF (max_abs > 1.0d-10) THEN
             rvec_e = MATMUL(at, DBLE(irvec_e(:, ir_e))) * alat * bohr2ang
             re_len = SQRT(DOT_PRODUCT(rvec_e, rvec_e))
             WRITE(iunit, '(2F14.6, ES16.8)') re_len, rp_len, max_abs
          ENDIF
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'Decay data written to ', TRIM(fname)
    WRITE(stdout, '(5X,A,I6,A,I6)') 'nrr_e=', nrr_e, '  nrr_p=', nrr_p
  END SUBROUTINE write_edmatw_2d_v2


  SUBROUTINE write_edmatw_2d_bin(nbndsub, nrr, irvec, ndegen, edmatw_2d, prefix_in)
    !-----------------------------------------------------------------------
    ! Write full M(R_e, R_p) in binary format for edwread restart.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: iunit, ir
    CHARACTER(LEN=256) :: fname

    iunit = 82
    fname = TRIM(prefix_in) // '_edmatw_2d.bin'
    OPEN(iunit, FILE=TRIM(fname), FORM='unformatted')
    WRITE(iunit) nbndsub, nrr
    WRITE(iunit) ndegen(1:nrr)
    DO ir = 1, nrr
       WRITE(iunit) irvec(:, ir)
    ENDDO
    WRITE(iunit) edmatw_2d
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'M(R_e,R_p) binary written to ', TRIM(fname)
  END SUBROUTINE write_edmatw_2d_bin


  SUBROUTINE read_edmatw_2d_file(prefix_in, nbndsub, nrr, ndegen, irvec, edmatw_2d)
    !-----------------------------------------------------------------------
    ! Read full M(R_e, R_p) from binary file.
    ! Ionode reads, then broadcasts to all ranks (avoids filesystem races).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in
    INTEGER, INTENT(OUT) :: nbndsub, nrr
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:), irvec(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw_2d(:,:,:,:)

    INTEGER :: iunit, ir
    CHARACTER(LEN=256) :: fname

    nbndsub = 0
    nrr = 0

    IF (ionode) THEN
       iunit = 82
       fname = TRIM(prefix_in) // '_edmatw_2d.bin'
       OPEN(iunit, FILE=TRIM(fname), FORM='unformatted', STATUS='old')
       READ(iunit) nbndsub, nrr
    ENDIF
    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr, ionode_id, world_comm)

    ALLOCATE(ndegen(nrr), irvec(3, nrr))

    IF (ionode) THEN
       READ(iunit) ndegen(1:nrr)
       DO ir = 1, nrr
          READ(iunit) irvec(:, ir)
       ENDDO
    ENDIF
    CALL mp_bcast(ndegen, ionode_id, world_comm)
    CALL mp_bcast(irvec, ionode_id, world_comm)

    ALLOCATE(edmatw_2d(nbndsub, nbndsub, nrr, nrr))

    IF (ionode) THEN
       READ(iunit) edmatw_2d
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'M(R_e,R_p) binary read from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I6)') '  nbndsub=', nbndsub, '  nrr=', nrr
    ENDIF
    CALL mp_bcast(edmatw_2d, ionode_id, world_comm)
  END SUBROUTINE read_edmatw_2d_file


  SUBROUTINE ed_fine_interp_2d(nbndsub, nrr, irvec, ndegen, chw, edmatw_2d, &
                                nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B with full double-FT M(R_e, R_p).
    !
    ! For each (k_i, k_f) on the fine grid:
    !   M_W(k_i,k_f) = Σ_{R_e,R_p} exp(ik_i·R_e) exp(-ik_f·R_p) M(R_e,R_p)
    !   Diagonalize H(k_i) → U(k_i), H(k_f) → U(k_f)
    !   M_B = U†(k_i) · M_W · U(k_f)
    !
    ! For diagonal (k_i = k_f): outputs prefix_edmatf_2d.dat
    ! For off-diagonal from K: outputs prefix_edmat_scatter_2d.dat
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch_2d
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, iunit_diag, iunit_scatter
    REAL(dp) :: xk_i(3), xk_f(3), eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_ki(nrr), cfac_kf(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub), edmatf_b(nbndsub, nbndsub)
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname
    REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
    COMPLEX(dp), PARAMETER :: ci = (0.0_dp, 1.0_dp)

    ! K point for off-diagonal scattering
    REAL(dp) :: xk_K(3)
    COMPLEX(dp) :: cfac_K(nrr), evec_K(nbndsub, nbndsub)
    REAL(dp) :: eig_K(nbndsub)

    nktotf = nk1f * nk2f * nk3f
    xk_K = (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp/)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B (2D): Interpolation using M(R_e, R_p)'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Diagonalize at K point once
    CALL get_cfac(nrr, irvec, xk_K, cfac_K)
    CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_K, evec_K, chw, cfac_K)

    iunit_diag = 80
    iunit_scatter = 79
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmatf_2d.dat'
       OPEN(iunit_diag, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit_diag, '(A)') '# Interpolated M from full double-FT M(R_e,R_p) — diagonal k_i=k_f'
       WRITE(iunit_diag, '(A)') '# ik  kx  ky  kz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'

       fname = TRIM(prefix_in) // '_edmat_scatter_2d.dat'
       OPEN(iunit_scatter, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit_scatter, '(A)') '# Interpolated M from full double-FT — K -> k_f'
       WRITE(iunit_scatter, '(A,3F10.5)') '# k_initial (crystal) = ', xk_K
       WRITE(iunit_scatter, '(A)') '# ik  kfx  kfy  kfz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_f(1) = DBLE(i1) / DBLE(nk1f)
             xk_f(2) = DBLE(i2) / DBLE(nk2f)
             xk_f(3) = DBLE(i3) / DBLE(nk3f)

             CALL get_cfac(nrr, irvec, xk_f, cfac_kf)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_kf)

             ! === Diagonal: k_i = k_f ===
             cfac_ki = cfac_kf
             CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_ki, cfac_kf, edmatf_w)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec_f, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit_diag, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, ABS(edmatf_b(ibnd,jbnd))**2, &
                           REAL(edmatf_b(ibnd,jbnd)), AIMAG(edmatf_b(ibnd,jbnd))
                   ENDDO
                ENDDO
             ENDIF

             ! === Off-diagonal: k_i = K, k_f varies ===
             CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_K, cfac_kf, edmatf_w)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec_K, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit_scatter, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, ABS(edmatf_b(ibnd,jbnd))**2, &
                           REAL(edmatf_b(ibnd,jbnd)), AIMAG(edmatf_b(ibnd,jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit_diag)
       CLOSE(iunit_scatter)
       WRITE(stdout, '(5X,A,I8,A)') 'Interpolated M(2D) at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A)') 'Diagonal: ' // TRIM(prefix_in) // '_edmatf_2d.dat'
       WRITE(stdout, '(5X,A)') 'K->k_f:  ' // TRIM(prefix_in) // '_edmat_scatter_2d.dat'
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp_2d

  SUBROUTINE write_edmatw(nbndsub, nrr, irvec, ndegen, edmatw, prefix_in)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir, ii, jj, iunit
    CHARACTER(LEN=256) :: fname

    iunit = 87
    fname = TRIM(prefix_in) // '_edmatw.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# Electron-defect matrix elements in Wannier basis'
    WRITE(iunit, '(I6)') nbndsub
    WRITE(iunit, '(I6)') nrr

    ir = 0
    DO WHILE (ir < nrr)
       WRITE(iunit, '(15I5)') ndegen(ir+1:MIN(ir+15, nrr))
       ir = MIN(ir + 15, nrr)
    ENDDO

    DO ir = 1, nrr
       DO jj = 1, nbndsub
          DO ii = 1, nbndsub
             WRITE(iunit, '(3I5, 2I5, 2E20.12)') irvec(:, ir), ii, jj, &
                  REAL(edmatw(ii, jj, ir), dp), AIMAG(edmatw(ii, jj, ir))
          ENDDO
       ENDDO
    ENDDO
    CLOSE(iunit)
  END SUBROUTINE write_edmatw

  SUBROUTINE read_edmatw_file(prefix_in, nbndsub, nrr, ndegen, irvec, edmatw)
    !-----------------------------------------------------------------------
    ! Read M(R) from prefix_edmatw.dat (written by write_edmatw).
    ! Allows skipping Part A and going directly to Part B.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: prefix_in
    INTEGER, INTENT(OUT) :: nbndsub, nrr
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:), irvec(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw(:,:,:)

    INTEGER :: iunit, ir, ii, jj, r1, r2, r3, idx_i, idx_j, ios
    INTEGER :: ntot, ierr
    REAL(dp) :: re_part, im_part
    REAL(dp), ALLOCATABLE :: rbuf(:)
    CHARACTER(LEN=256) :: fname, line

    iunit = 88
    fname = TRIM(prefix_in) // '_edmatw.dat'
    nbndsub = 0; nrr = 0

    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('read_edmatw_file', 'Cannot open ' // TRIM(fname), 1)

       READ(iunit, '(A)') line  ! header
       READ(iunit, *) nbndsub
       READ(iunit, *) nrr

       ALLOCATE(ndegen(nrr))
       ir = 0
       DO WHILE (ir < nrr)
          READ(iunit, *, IOSTAT=ios) ndegen(ir+1:MIN(ir+15, nrr))
          IF (ios /= 0) EXIT
          ir = MIN(ir + 15, nrr)
       ENDDO

       ALLOCATE(irvec(3, nrr))
       ALLOCATE(edmatw(nbndsub, nbndsub, nrr))
       edmatw = (0.0_dp, 0.0_dp)

       DO ir = 1, nrr
          DO jj = 1, nbndsub
             DO ii = 1, nbndsub
                READ(iunit, *) r1, r2, r3, idx_i, idx_j, re_part, im_part
                IF (ii == 1 .AND. jj == 1) irvec(:, ir) = (/r1, r2, r3/)
                edmatw(idx_i, idx_j, ir) = CMPLX(re_part, im_part, dp)
             ENDDO
          ENDDO
       ENDDO

       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Read M(R) from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I6)') 'nbndsub = ', nbndsub, ', nrr = ', nrr
    ENDIF

    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr, ionode_id, world_comm)

    IF (.NOT. ionode) THEN
       ALLOCATE(ndegen(nrr))
       ALLOCATE(irvec(3, nrr))
       ALLOCATE(edmatw(nbndsub, nbndsub, nrr))
    ENDIF

    CALL mp_bcast(ndegen, ionode_id, world_comm)
    CALL mp_bcast(irvec, ionode_id, world_comm)

    ! Broadcast complex 3D array as two real 1D arrays
    ntot = nbndsub * nbndsub * nrr
    ALLOCATE(rbuf(ntot))
    IF (ionode) rbuf = RESHAPE(REAL(edmatw, dp), (/ntot/))
    CALL MPI_Bcast(rbuf, ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) edmatw = RESHAPE(CMPLX(rbuf, 0.0_dp, dp), (/nbndsub, nbndsub, nrr/))
    IF (ionode) rbuf = RESHAPE(AIMAG(edmatw), (/ntot/))
    CALL MPI_Bcast(rbuf, ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) &
         edmatw = CMPLX(REAL(edmatw, dp), RESHAPE(rbuf, (/nbndsub, nbndsub, nrr/)), dp)
    DEALLOCATE(rbuf)

  END SUBROUTINE read_edmatw_file

  SUBROUTINE load_supercell_pot(prefix_in, outdir_in, vf_out)
    !-----------------------------------------------------------------------
    ! Load supercell potential: ionode reads serially, broadcasts to all.
    !
    ! Strategy:
    !   1. ALL ranks collectively clean_pw(.TRUE.) to free memory
    !   2. ionode sets communicators to MPI_COMM_SELF and calls
    !      get_vloc_onthefly (serial: only ~2GB needed on 1 rank)
    !   3. ionode cleans up supercell QE state
    !   4. ionode broadcasts the extracted potential to all ranks
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp_world, ONLY : world_comm
    USE mp_pools, ONLY : intra_pool_comm, inter_pool_comm, npool, my_pool_id
    USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm, nbgrp
    USE mp_images, ONLY : intra_image_comm, nimage
    USE mp, ONLY : mp_bcast, mp_barrier
    USE io_files, ONLY : prefix, tmp_dir
    USE edic_mod, ONLY : V_file
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: prefix_in, outdir_in
    TYPE(V_file), INTENT(INOUT) :: vf_out

    INTERFACE
       SUBROUTINE get_vloc_onthefly(prefix_in, outdir_in, vf_out, is_noncolin, bxc_out)
         USE kinds, ONLY : dp
         USE edic_mod, ONLY : V_file
         CHARACTER(LEN=*), INTENT(IN) :: prefix_in, outdir_in
         TYPE(V_file), INTENT(INOUT) :: vf_out
         LOGICAL, INTENT(IN) :: is_noncolin
         REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: bxc_out(:,:)
       END SUBROUTINE
    END INTERFACE

    INTEGER :: save_world, save_pool, save_inter_pool
    INTEGER :: save_bgrp, save_inter_bgrp, save_image
    INTEGER :: save_npool, save_pool_id, save_nbgrp, save_nimage
    INTEGER :: nrtot, ierr
    REAL(dp), ALLOCATABLE :: pot_tmp(:)
    CHARACTER(LEN=256) :: save_prefix, save_tmpdir

    ! Save ALL state
    save_world = world_comm
    save_pool = intra_pool_comm
    save_inter_pool = inter_pool_comm
    save_bgrp = intra_bgrp_comm
    save_inter_bgrp = inter_bgrp_comm
    save_image = intra_image_comm
    save_npool = npool
    save_pool_id = my_pool_id
    save_nbgrp = nbgrp
    save_nimage = nimage
    save_prefix = prefix
    save_tmpdir = tmp_dir

    ! Step 1: ALL ranks collectively clean QE state (free memory)
    CALL clean_pw(.TRUE.)

    ! Step 2: ionode reads supercell potential serially
    IF (ionode) THEN
       world_comm = MPI_COMM_SELF
       intra_pool_comm = MPI_COMM_SELF
       inter_pool_comm = MPI_COMM_SELF
       intra_bgrp_comm = MPI_COMM_SELF
       inter_bgrp_comm = MPI_COMM_SELF
       intra_image_comm = MPI_COMM_SELF
       npool = 1
       my_pool_id = 0
       nbgrp = 1
       nimage = 1

       CALL get_vloc_onthefly(prefix_in, outdir_in, vf_out, .FALSE.)
       CALL clean_pw(.TRUE.)
    ENDIF

    ! Step 3: Restore ALL state on ALL ranks
    world_comm = save_world
    intra_pool_comm = save_pool
    inter_pool_comm = save_inter_pool
    intra_bgrp_comm = save_bgrp
    inter_bgrp_comm = save_inter_bgrp
    intra_image_comm = save_image
    npool = save_npool
    my_pool_id = save_pool_id
    nbgrp = save_nbgrp
    nimage = save_nimage
    prefix = save_prefix
    tmp_dir = save_tmpdir

    ! Step 4: Broadcast potential data from ionode to all ranks
    ! Use MPI_Bcast with explicit counts — ionode's pot array may be
    ! larger than nr1*nr2*nr3 due to FFT padding (dfftp%nnr > nr1*nr2*nr3)
    CALL MPI_Bcast(vf_out%nr1, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nr2, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nr3, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nr1x, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nr2x, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nr3x, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%nat, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%ntyp, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%alat, 1, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%ibrav, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%celldm, 6, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%at, 9, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%omega, 1, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)

    nrtot = vf_out%nr1 * vf_out%nr2 * vf_out%nr3

    ! On ionode, pot may have FFT padding (size > nrtot). Trim to nrtot.
    IF (ionode) THEN
       ALLOCATE(pot_tmp(nrtot))
       pot_tmp(:) = vf_out%pot(1:nrtot)
       DEALLOCATE(vf_out%pot)
       ALLOCATE(vf_out%pot(nrtot))
       vf_out%pot(:) = pot_tmp(:)
       DEALLOCATE(pot_tmp)
    ELSE
       ALLOCATE(vf_out%pot(nrtot))
       ALLOCATE(vf_out%ityp(vf_out%nat))
       ALLOCATE(vf_out%tau(3, vf_out%nat))
    ENDIF
    CALL MPI_Bcast(vf_out%pot, nrtot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%ityp, vf_out%nat, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(vf_out%tau, 3*vf_out%nat, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)

    IF (ionode) WRITE(stdout, '(5X,A,3I6,A,I6)') &
         '  grid: ', vf_out%nr1, vf_out%nr2, vf_out%nr3, &
         '  nat:', vf_out%nat

  END SUBROUTINE load_supercell_pot

  SUBROUTINE read_filukk_edi(fname, nbnd, nkstot, nks, nbndsub)
    !-----------------------------------------------------------------------
    ! Read Wannier rotation matrices from filukk file.
    ! Sets u_mat, u_mat_opt, lwindow, excluded_band, num_bands, n_wannier
    ! in wann_common module, so ed_coarse_calc can use them.
    ! The filukk stores the COMBINED rotation u_kc = u_opt * u.
    ! We store u_kc in u_mat and set num_bands = n_wannier so
    ! ed_coarse_calc uses cu = u_mat directly (no u_mat_opt needed).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    USE wann_common, ONLY : u_mat, u_mat_opt, lwindow, excluded_band, &
                             n_wannier, num_bands, iknum, wann_centers
    USE global_var, ONLY : nbndep, nbndskip, ibndkept
    USE parallelism, ONLY : fkbounds
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: nbnd, nkstot, nks, nbndsub

    INTEGER :: ik, ibnd, iw, iunit, nbndep_f, nbndskip_f
    INTEGER :: lower_bnd, upper_bnd, ntot, ierr
    COMPLEX(dp), ALLOCATABLE :: u_kc_all(:,:,:)
    REAL(dp), ALLOCATABLE :: rbuf(:)
    LOGICAL, ALLOCATABLE :: lwindow_all(:,:)
    INTEGER, ALLOCATABLE :: lwin_int(:)

    iunit = 89
    iknum = nkstot

    ! nbndsub may not be broadcast to non-ionode ranks yet.
    ! Use ionode's value and broadcast.
    IF (ionode) n_wannier = nbndsub
    CALL MPI_Bcast(n_wannier, 1, MPI_INTEGER, ionode_id, world_comm, ierr)

    IF (ionode) THEN
       OPEN(UNIT=iunit, FILE=TRIM(fname), FORM='formatted', STATUS='old')

       ! Read nbndep, nbndskip
       READ(iunit, *) nbndep_f, nbndskip_f
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(nbndep_f, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(nbndskip_f, 1, MPI_INTEGER, ionode_id, world_comm, ierr)

    nbndep = nbndep_f
    nbndskip = nbndskip_f
    num_bands = nbndep

    ! Read ibndkept
    IF (.NOT. ALLOCATED(ibndkept)) ALLOCATE(ibndkept(nbndep))
    IF (ionode) THEN
       DO ibnd = 1, nbndep
          READ(iunit, *) ibndkept(ibnd)
       ENDDO
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(ibndkept, nbndep, MPI_INTEGER, ionode_id, world_comm, ierr)

    ! Read u_kc (combined rotation) for all k-points
    ALLOCATE(u_kc_all(nbndep, n_wannier, nkstot))
    u_kc_all = (0.0_dp, 0.0_dp)
    IF (ionode) THEN
       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             DO iw = 1, n_wannier
                READ(iunit, *) u_kc_all(ibnd, iw, ik)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! Broadcast u_kc using direct MPI_Bcast (QE wrapper has issues with large arrays)
    ntot = nbndep * n_wannier * nkstot
    ALLOCATE(rbuf(2*ntot))
    IF (ionode) THEN
       DO iw = 1, ntot
          rbuf(2*iw-1) = REAL(u_kc_all(MOD(iw-1,nbndep)+1, &
               MOD((iw-1)/nbndep,n_wannier)+1, (iw-1)/(nbndep*n_wannier)+1), dp)
          rbuf(2*iw) = AIMAG(u_kc_all(MOD(iw-1,nbndep)+1, &
               MOD((iw-1)/nbndep,n_wannier)+1, (iw-1)/(nbndep*n_wannier)+1))
       ENDDO
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(rbuf, 2*ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) THEN
       DO iw = 1, ntot
          u_kc_all(MOD(iw-1,nbndep)+1, MOD((iw-1)/nbndep,n_wannier)+1, &
               (iw-1)/(nbndep*n_wannier)+1) = CMPLX(rbuf(2*iw-1), rbuf(2*iw), dp)
       ENDDO
    ENDIF
    DEALLOCATE(rbuf)

    ! Read lwindow for all k-points (broadcast as integer via MPI_Bcast)
    ALLOCATE(lwindow_all(nbndep, nkstot))
    ntot = nbndep * nkstot
    ALLOCATE(lwin_int(ntot))
    lwin_int = 0
    IF (ionode) THEN
       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             READ(iunit, *) lwindow_all(ibnd, ik)
             IF (lwindow_all(ibnd, ik)) lwin_int((ik-1)*nbndep + ibnd) = 1
          ENDDO
       ENDDO
    ENDIF
    CALL MPI_Bcast(lwin_int, ntot, MPI_INTEGER, ionode_id, world_comm, ierr)
    DO ik = 1, nkstot
       DO ibnd = 1, nbndep
          lwindow_all(ibnd, ik) = lwin_int((ik-1)*nbndep + ibnd) > 0
       ENDDO
    ENDDO
    DEALLOCATE(lwin_int)

    ! Read excluded_band (broadcast as integer via MPI_Bcast)
    IF (.NOT. ALLOCATED(excluded_band)) ALLOCATE(excluded_band(nbnd))
    ALLOCATE(lwin_int(nbnd))
    lwin_int = 0
    IF (ionode) THEN
       DO ibnd = 1, nbnd
          READ(iunit, *) excluded_band(ibnd)
          IF (excluded_band(ibnd)) lwin_int(ibnd) = 1
       ENDDO
    ENDIF
    CALL MPI_Bcast(lwin_int, nbnd, MPI_INTEGER, ionode_id, world_comm, ierr)
    excluded_band = lwin_int > 0
    DEALLOCATE(lwin_int)

    ! Read wann_centers
    IF (.NOT. ALLOCATED(wann_centers)) ALLOCATE(wann_centers(3, n_wannier))
    IF (ionode) THEN
       DO iw = 1, n_wannier
          READ(iunit, *) wann_centers(:, iw)
       ENDDO
       CLOSE(iunit)
    ENDIF
    CALL MPI_Bcast(wann_centers, 3*n_wannier, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)

    ! Extract pool-local k-points
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)

    ! Store u_kc in u_mat (pool-local), set num_bands = n_wannier
    ! so ed_coarse_calc uses cu = u_mat directly
    IF (ALLOCATED(u_mat)) DEALLOCATE(u_mat)
    ALLOCATE(u_mat(nbndep, n_wannier, nks))
    u_mat(:,:,1:nks) = u_kc_all(:,:,lower_bnd:upper_bnd)

    ! u_mat_opt not needed (num_bands = n_wannier path), but allocate for safety
    IF (ALLOCATED(u_mat_opt)) DEALLOCATE(u_mat_opt)
    ALLOCATE(u_mat_opt(nbndep, nbndep, nks))
    u_mat_opt = (0.0_dp, 0.0_dp)

    ! lwindow (pool-local)
    IF (ALLOCATED(lwindow)) DEALLOCATE(lwindow)
    ALLOCATE(lwindow(nbndep, nks))
    lwindow(:,1:nks) = lwindow_all(:,lower_bnd:upper_bnd)

    DEALLOCATE(u_kc_all, lwindow_all)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,A)') 'Read Wannier data from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I4,A,I4)') &
            'nbndep=', nbndep, ' n_wannier=', n_wannier, ' nkstot=', nkstot
    ENDIF

  END SUBROUTINE read_filukk_edi


  SUBROUTINE generate_nscf_input(filki, filkf, edi_outdir_in)
    !-----------------------------------------------------------------------
    ! Auto-generate nscf_custom.in for k-points not in the current NSCF data.
    ! Reads parameters from QE modules (populated by read_file).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE io_files, ONLY : prefix, pseudo_dir
    USE cell_base, ONLY : at, alat, tpiba2
    USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, atm
    USE wvfct, ONLY : nbnd
    USE gvect, ONLY : gcutm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, edi_outdir_in
    REAL(dp) :: ecutwfc_val, ecutrho_val

    INTEGER :: nki, nkf, nk_merged, ik, jk, iunit
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:), xk_merged(:,:)
    REAL(dp) :: diff(3), dist
    LOGICAL :: is_dup
    INTEGER :: ia

    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (.NOT. ionode) THEN
       IF (ALLOCATED(xki)) DEALLOCATE(xki)
       IF (ALLOCATED(xkf)) DEALLOCATE(xkf)
       RETURN
    ENDIF

    ! Merge ki and kf lists, removing duplicates
    ALLOCATE(xk_merged(3, nki + nkf))
    nk_merged = 0
    DO ik = 1, nki
       is_dup = .FALSE.
       DO jk = 1, nk_merged
          diff = xki(:,ik) - xk_merged(:,jk)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-10) THEN
             is_dup = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF (.NOT. is_dup) THEN
          nk_merged = nk_merged + 1
          xk_merged(:, nk_merged) = xki(:, ik)
       ENDIF
    ENDDO
    DO ik = 1, nkf
       is_dup = .FALSE.
       DO jk = 1, nk_merged
          diff = xkf(:,ik) - xk_merged(:,jk)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-10) THEN
             is_dup = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF (.NOT. is_dup) THEN
          nk_merged = nk_merged + 1
          xk_merged(:, nk_merged) = xkf(:, ik)
       ENDIF
    ENDDO

    ! Compute cutoffs from QE internal variables (Ry units)
    ! In QE: ecutwfc = ecfixed or from input; gcutm = ecutrho/tpiba2
    ! We recover ecutrho = gcutm * tpiba2
    ecutrho_val = gcutm * tpiba2
    ! ecutwfc: use the ratio ecutrho/ecutwfc which is typically 4 (NC) or 8-12 (US/PAW)
    ! We can get ecutwfc from gvecw module
    BLOCK
       USE gvecw, ONLY : ecutwfc_in => ecutwfc
       ecutwfc_val = ecutwfc_in
    END BLOCK

    ! Write nscf_custom.in
    iunit = 90
    OPEN(iunit, FILE='nscf_custom.in', FORM='formatted')
    WRITE(iunit, '(A)')       '&CONTROL'
    WRITE(iunit, '(A)')       "  calculation = 'bands'"
    WRITE(iunit, '(A,A,A)')   "  prefix = '", TRIM(prefix), "'"
    WRITE(iunit, '(A,A,A)')   "  outdir = '", TRIM(edi_outdir_in), "'"
    WRITE(iunit, '(A,A,A)')   "  pseudo_dir = '", TRIM(pseudo_dir), "'"
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       '&SYSTEM'
    WRITE(iunit, '(A)')       '  ibrav = 0'
    WRITE(iunit, '(A,F12.6)') '  ecutwfc = ', ecutwfc_val
    WRITE(iunit, '(A,F12.6)') '  ecutrho = ', ecutrho_val
    WRITE(iunit, '(A,I4)')    '  nat = ', nat
    WRITE(iunit, '(A,I4)')    '  ntyp = ', ntyp
    WRITE(iunit, '(A,I4)')    '  nbnd = ', nbnd
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       '&ELECTRONS'
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       'ATOMIC_SPECIES'
    DO ia = 1, ntyp
       WRITE(iunit, '(A,A,A)') TRIM(atm(ia)), '  0.0  ', TRIM(atm(ia)) // '.upf'
    ENDDO
    WRITE(iunit, '(A)')       'CELL_PARAMETERS alat'
    WRITE(iunit, '(3F16.10)') at(:,1)
    WRITE(iunit, '(3F16.10)') at(:,2)
    WRITE(iunit, '(3F16.10)') at(:,3)
    WRITE(iunit, '(A)')       'ATOMIC_POSITIONS crystal'
    DO ia = 1, nat
       WRITE(iunit, '(A,3F16.10)') TRIM(atm(ityp(ia))), tau(:,ia)
    ENDDO
    WRITE(iunit, '(A)')       'K_POINTS crystal'
    WRITE(iunit, '(I6)')      nk_merged
    DO ik = 1, nk_merged
       WRITE(iunit, '(3F16.10,A)') xk_merged(:,ik), '  1.0'
    ENDDO
    CLOSE(iunit)

    WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
    WRITE(stdout, '(5X,A)')   'NSCF wavefunctions not found for required k-points.'
    WRITE(stdout, '(5X,A,I6,A)') 'Auto-generated: nscf_custom.in (', nk_merged, ' k-points)'
    WRITE(stdout, '(5X,A)')   ''
    WRITE(stdout, '(5X,A)')   'Please run:  pw.x < nscf_custom.in > nscf_custom.out'
    WRITE(stdout, '(5X,A)')   'Then re-run EDI with the same input.'
    WRITE(stdout, '(5X,A)')   REPEAT('=', 60)

    DEALLOCATE(xki, xkf, xk_merged)
  END SUBROUTINE generate_nscf_input


  SUBROUTINE ed_direct_from_files(filki, filkf, prefix_in, band_ed_str)
    !-----------------------------------------------------------------------
    ! Direct calculation mode: compute M(k_i, k_f) from wavefunctions
    ! for arbitrary k-points read from files.
    !
    ! band_ed_str: band selection string, e.g. '13-17' for bands 13 to 17.
    !   If empty, all bands are used.
    !
    ! Prerequisites: primitive cell state restored (read_file_new done),
    !   supercell potentials loaded in V_d, V_p, V_colin.
    ! Requires NSCF wavefunctions at ALL requested k-points.
    ! Supercell potentials are loaded internally.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast, mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_global, ONLY : inter_pool_comm
    USE mp_pools, ONLY : npool, my_pool_id
    USE cell_base, ONLY : at
    USE klist, ONLY : xk, igk_k, ngk, nkstot, nks
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module, ONLY : npol
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE edic_mod, ONLY : V_d, V_colin
    USE ep_constants, ONLY : czero, cone
    USE constants, ONLY : tpi
    USE parallelism, ONLY : fkbounds
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, prefix_in, band_ed_str

    INTEGER :: nki, nkf, iki, ikf, ik_nscf, ibnd, jbnd, ierr_mpi
    INTEGER :: ig, ir, inr, irx, iry, irz, ir1mod, ir2mod, ir3mod, irnmod
    INTEGER :: iunit, lower_bnd, upper_bnd, nbnd_kept
    INTEGER :: ibnd_min, ibnd_max, ipos, ios_parse
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:), xk_cryst_nscf(:,:)
    REAL(dp) :: qcryst(3), d1, d2, d3, arg, diff(3)
    COMPLEX(dp) :: phase, mlocal
    COMPLEX(dp), ALLOCATABLE :: psir_ki(:,:,:), psir_kf(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: evc_tmp(:,:), psic_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: V_folded(:), edmatkq(:,:)
    INTEGER :: ipol_d
    INTEGER, ALLOCATABLE :: iki_map(:), ikf_map(:)
    CHARACTER(LEN=256) :: fname

    ! Read k-point files
    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Direct calculation: M(k_i, k_f) from wavefunctions'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6,A,I6,A,I10)') &
            'nki=', nki, '  nkf=', nkf, '  pairs=', nki * nkf
       FLUSH(stdout)
    ENDIF

    ! Build crystal-coord map of NSCF k-points
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)
    ALLOCATE(xk_cryst_nscf(3, nkstot))
    xk_cryst_nscf = 0.0_dp
    xk_cryst_nscf(:, lower_bnd:upper_bnd) = xk(:, 1:nks)
    CALL cryst_to_cart(nks, xk_cryst_nscf(:, lower_bnd:upper_bnd), at, -1)
    CALL mp_sum(xk_cryst_nscf, inter_pool_comm)

    ! Map user k-points to NSCF indices
    ALLOCATE(iki_map(nki), ikf_map(nkf))
    iki_map = 0
    ikf_map = 0

    DO iki = 1, nki
       DO ik_nscf = 1, nkstot
          diff = xki(:,iki) - xk_cryst_nscf(:,ik_nscf)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-8) THEN
             iki_map(iki) = ik_nscf
             EXIT
          ENDIF
       ENDDO
    ENDDO
    DO ikf = 1, nkf
       DO ik_nscf = 1, nkstot
          diff = xkf(:,ikf) - xk_cryst_nscf(:,ik_nscf)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-8) THEN
             ikf_map(ikf) = ik_nscf
             EXIT
          ENDIF
       ENDDO
    ENDDO

    ! Check if all k-points are found
    IF (ANY(iki_map == 0) .OR. ANY(ikf_map == 0)) THEN
       IF (ionode) THEN
          WRITE(stdout, '(5X,A)') 'Some k-points not found in NSCF data.'
          DO iki = 1, nki
             IF (iki_map(iki) == 0) WRITE(stdout, '(5X,A,3F10.5)') &
                  '  Missing ki: ', xki(:,iki)
          ENDDO
          DO ikf = 1, nkf
             IF (ikf_map(ikf) == 0) WRITE(stdout, '(5X,A,3F10.5)') &
                  '  Missing kf: ', xkf(:,ikf)
          ENDDO
       ENDIF
       ! Auto-generate NSCF input and stop
       BLOCK
          USE io_files, ONLY : tmp_dir_loc => tmp_dir
          CALL generate_nscf_input(filki, filkf, tmp_dir_loc)
       END BLOCK
       DEALLOCATE(xki, xkf, xk_cryst_nscf, iki_map, ikf_map)
       CALL errore('ed_direct_from_files', &
            'k-points not in NSCF data. Run pw.x with nscf_custom.in first.', 1)
    ENDIF

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'All k-points found in NSCF data.'
       FLUSH(stdout)
    ENDIF

    ! Parse band_ed string (e.g. '13-17') to get band range
    IF (LEN_TRIM(band_ed_str) > 0) THEN
       ipos = INDEX(band_ed_str, '-')
       IF (ipos > 0) THEN
          READ(band_ed_str(1:ipos-1), *, IOSTAT=ios_parse) ibnd_min
          READ(band_ed_str(ipos+1:), *, IOSTAT=ios_parse) ibnd_max
       ELSE
          READ(band_ed_str, *, IOSTAT=ios_parse) ibnd_min
          ibnd_max = ibnd_min
       ENDIF
       IF (ibnd_min < 1) ibnd_min = 1
       IF (ibnd_max > nbnd) ibnd_max = nbnd
       IF (ionode) WRITE(stdout, '(5X,A,I4,A,I4,A,I4,A)') &
            'band_ed = ', ibnd_min, ' - ', ibnd_max, &
            '  (', ibnd_max - ibnd_min + 1, ' bands selected)'
    ELSE
       ibnd_min = 1
       ibnd_max = nbnd
       IF (ionode) WRITE(stdout, '(5X,A,I4,A)') &
            'band_ed not set, using all ', nbnd, ' bands'
    ENDIF
    nbnd_kept = ibnd_max - ibnd_min + 1

    ! Potentials (V_colin) and primitive cell state already set up by caller
    ! Panel broadcast: cache local wfcs, broadcast per pool (same as double-FT)
    ALLOCATE(psic_tmp(dffts%nnr))
    ALLOCATE(evc_tmp(npwx * npol, nbnd))
    ALLOCATE(V_folded(dffts%nnr))
    ALLOCATE(edmatkq(nbnd_kept, nbnd_kept))
    ALLOCATE(psir_ki(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psir_kf(dffts%nnr, npol, nbnd_kept))

    BLOCK
       COMPLEX(dp), ALLOCATABLE :: psir_cache(:,:,:,:), psir_recv(:,:,:,:)
       INTEGER :: ik_cache, ib_cache, ig_cache, npairs_done, ipol_c
       INTEGER :: ip, src_lower, src_nks, nkbase, nkrest, nks_max
       INTEGER :: ikf_local, ik_pool, ik_l

       nkbase = nkstot / npool
       nkrest = MOD(nkstot, npool)
       nks_max = nkbase
       IF (nkrest > 0) nks_max = nkbase + 1

       ! Step 1: Cache pool-local wavefunctions in real space (read disk ONCE)
       ! Only bands ibnd_min..ibnd_max are stored (index 1..nbnd_kept)
       ALLOCATE(psir_cache(dffts%nnr, npol, nbnd_kept, nks))
       DO ik_cache = 1, nks
          CALL read_collected_wfc(restart_dir(), ik_cache, evc_tmp)
          DO ib_cache = 1, nbnd_kept
             DO ipol_c = 1, npol
                psic_tmp = (0.0_dp, 0.0_dp)
                DO ig_cache = 1, ngk(ik_cache)
                   psic_tmp(dffts%nl(igk_k(ig_cache, ik_cache))) = &
                        evc_tmp(ig_cache + (ipol_c-1)*npwx, ibnd_min + ib_cache - 1)
                ENDDO
                CALL invfft('Wave', psic_tmp, dffts)
                psir_cache(:, ipol_c, ib_cache, ik_cache) = psic_tmp
             ENDDO
          ENDDO
       ENDDO

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,I4,A)') &
               'Cached ', nks, ' pool-local k-points in real space'
          FLUSH(stdout)
       ENDIF

       ! Step 2: Get psi(ki) — find in local cache or receive via panel broadcast
       ! For now, broadcast ki from owning pool (only nki broadcasts, typically 1)
       DO iki = 1, nki
          CALL get_pool_and_local(iki_map(iki), nkstot, npool, ik_pool, ik_l)
          psir_ki = (0.0_dp, 0.0_dp)
          IF (my_pool_id == ik_pool) psir_ki(:,:,:) = psir_cache(:,:,:,ik_l)
          CALL MPI_Bcast(psir_ki, dffts%nnr * npol * nbnd_kept * 2, &
               MPI_DOUBLE_PRECISION, ik_pool, inter_pool_comm, ierr_mpi)
       ENDDO
       ! (psir_ki now holds the last ki — fine for nki=1)

       ! Step 3: Panel broadcast for kf — iterate over source pools
       ALLOCATE(psir_recv(dffts%nnr, npol, nbnd_kept, nks_max))

       iunit = 87
       IF (ionode) THEN
          fname = TRIM(prefix_in) // '_edmat_direct.dat'
          OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
          WRITE(iunit, '(A)') '# Direct e-d matrix elements M(n,m, k_i, k_f)'
          WRITE(iunit, '(A,I6,A,I6)') '# nki=', nki, ' nkf=', nkf
          WRITE(iunit, '(A,I4,A,I4)') '# band range: ', ibnd_min, ' - ', ibnd_max
          WRITE(iunit, '(A)') '# iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)'
       ENDIF

       npairs_done = 0

       DO ip = 0, npool - 1
          ! Determine k-point range of source pool ip
          IF (ip < nkrest) THEN
             src_lower = ip * (nkbase + 1) + 1
             src_nks = nkbase + 1
          ELSE
             src_lower = nkrest * (nkbase + 1) + (ip - nkrest) * nkbase + 1
             src_nks = nkbase
          ENDIF

          ! Source pool broadcasts its cached wavefunctions
          IF (my_pool_id == ip) THEN
             psir_recv(:, :, :, 1:src_nks) = psir_cache(:, :, :, 1:nks)
          ENDIF
          CALL MPI_Bcast(psir_recv, dffts%nnr * npol * nbnd_kept * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)

          ! Process all kf that map to source pool ip's k-points
          DO iki = 1, nki
             DO ikf = 1, nkf
                ! Check if this kf's NSCF index is in source pool ip's range
                IF (ikf_map(ikf) < src_lower .OR. &
                    ikf_map(ikf) >= src_lower + src_nks) CYCLE

                ikf_local = ikf_map(ikf) - src_lower + 1
                psir_kf(:,:,:) = psir_recv(:,:,:, ikf_local)

                ! Compute V_folded and matrix element
                qcryst = xkf(:,ikf) - xki(:,iki)
                V_folded = (0.0_dp, 0.0_dp)
                IF (ABS(qcryst(1)) < 1.0d-8 .AND. ABS(qcryst(2)) < 1.0d-8 &
                     .AND. ABS(qcryst(3)) < 1.0d-8) THEN
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
                            V_folded(irnmod) = V_folded(irnmod) + V_colin(inr)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   d1 = tpi * qcryst(1) / dffts%nr1
                   d2 = tpi * qcryst(2) / dffts%nr2
                   d3 = tpi * qcryst(3) / dffts%nr3
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
                            arg = irx * d1 + iry * d2 + irz * d3
                            phase = CMPLX(COS(arg), SIN(arg), dp)
                            V_folded(irnmod) = V_folded(irnmod) + V_colin(inr) * phase
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

                ! Local matrix element: sum over spinor components
                DO ibnd = 1, nbnd_kept
                   DO jbnd = 1, nbnd_kept
                      mlocal = czero
                      DO ipol_d = 1, npol
                         DO ir = 1, dffts%nnr
                            mlocal = mlocal + CONJG(psir_ki(ir, ipol_d, ibnd)) * &
                                 psir_kf(ir, ipol_d, jbnd) * V_folded(ir)
                         ENDDO
                      ENDDO
                      edmatkq(ibnd, jbnd) = mlocal / DBLE(dffts%nnr)
                   ENDDO
                ENDDO

                npairs_done = npairs_done + 1
                IF (ionode) THEN
                   DO ibnd = 1, nbnd_kept
                      DO jbnd = 1, nbnd_kept
                         WRITE(iunit, '(2I6,6F10.5,2I4,3ES16.8)') &
                              iki, ikf, xki(:,iki), xkf(:,ikf), &
                              ibnd_min + ibnd - 1, ibnd_min + jbnd - 1, &
                              ABS(edmatkq(ibnd, jbnd))**2, &
                              REAL(edmatkq(ibnd, jbnd)), AIMAG(edmatkq(ibnd, jbnd))
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO  ! ikf
          ENDDO  ! iki

          IF (ionode) THEN
             WRITE(stdout, '(5X,A,I4,A,I4,A,I6,A,I6)') &
                  'kf pool ', ip, ' / ', npool - 1, &
                  '  pairs: ', npairs_done, ' / ', nki * nkf
             FLUSH(stdout)
          ENDIF
       ENDDO  ! ip (source pool)

       IF (ionode) THEN
          CLOSE(iunit)
          fname = TRIM(prefix_in) // '_edmat_direct.dat'
          WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
          WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       ENDIF

       DEALLOCATE(psir_cache, psir_recv)
    END BLOCK

    DEALLOCATE(xki, xkf, xk_cryst_nscf, iki_map, ikf_map)
    DEALLOCATE(evc_tmp, psic_tmp, psir_ki, psir_kf, V_folded, edmatkq)

  CONTAINS
    SUBROUTINE get_pool_and_local(ik_g, nktot, np, ipool, ik_l)
      ! Map global k-index to owning pool and local index
      INTEGER, INTENT(IN) :: ik_g, nktot, np
      INTEGER, INTENT(OUT) :: ipool, ik_l
      INTEGER :: nkb, nkr, off
      nkb = nktot / np
      nkr = MOD(nktot, np)
      IF (ik_g <= (nkb + 1) * nkr) THEN
         ipool = (ik_g - 1) / (nkb + 1)
         ik_l = ik_g - ipool * (nkb + 1)
      ELSE
         ipool = nkr + (ik_g - (nkb + 1) * nkr - 1) / nkb
         off = nkr * (nkb + 1) + (ipool - nkr) * nkb
         ik_l = ik_g - off
      ENDIF
    END SUBROUTINE
  END SUBROUTINE ed_direct_from_files


  SUBROUTINE read_kpoints_file(fname, nk, xk_cryst)
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(OUT) :: nk
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: xk_cryst(:,:)
    INTEGER :: ik, iunit, ios

    iunit = 84
    nk = 0
    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('read_kpoints_file', 'Cannot open ' // TRIM(fname), 1)
       READ(iunit, *) nk
    ENDIF
    CALL mp_bcast(nk, ionode_id, world_comm)
    ALLOCATE(xk_cryst(3, nk))
    IF (ionode) THEN
       DO ik = 1, nk
          READ(iunit, *) xk_cryst(:, ik)
       ENDDO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A,A,I6,A)') 'Read ', TRIM(fname), ': ', nk, ' k-points'
    ENDIF
    CALL mp_bcast(xk_cryst, ionode_id, world_comm)
  END SUBROUTINE read_kpoints_file


  SUBROUTINE ed_interp_from_file(nbndsub, nrr, irvec, ndegen, chw, edmatw_2d, &
                                  filki, filkf, prefix_in)
    !-----------------------------------------------------------------------
    ! Wannier interpolation from k-point files using full double-FT
    ! M(R_e, R_p) via edmatwan2bloch_2d.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch_2d
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, prefix_in

    INTEGER :: iki, ikf, nki, nkf, ibnd, jbnd, iunit
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:)
    REAL(dp) :: eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_ki(nrr), cfac_kf(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub), tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname

    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)') 'Wannier interpolation from k-point files [double-FT]'
       WRITE(stdout, '(5X,A)') '  Using M(R_e, R_p) with edmatwan2bloch_2d'
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6,A,I6,A,I10)') &
            'nki=', nki, ' nkf=', nkf, ' pairs=', nki*nkf

       iunit = 83
       fname = TRIM(prefix_in) // '_edmat_interp.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Wannier-interpolated M(n,m, k_i, k_f) [double-FT]'
       WRITE(iunit, '(A,I6,A,I6)') '# nki=', nki, ' nkf=', nkf
       WRITE(iunit, '(A)') '# iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    DO iki = 1, nki
       ! Diag H(ki) -> U(ki)
       CALL get_cfac(nrr, irvec, xki(:,iki), cfac_ki)
       CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_i, evec_i, chw, cfac_ki)

       DO ikf = 1, nkf
          ! Diag H(kf) -> U(kf)
          CALL get_cfac(nrr, irvec, xkf(:,ikf), cfac_kf)
          CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_kf)

          ! Double-FT inverse: M_W(ki,kf) from M(Re, Rp)
          CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_ki, cfac_kf, edmatf_w)

          ! Rotate to Bloch: M_B = U_dag(ki) M_W U(kf)
          CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                      cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
          CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                      cone, evec_i, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

          IF (ionode) THEN
             DO ibnd = 1, nbndsub
                DO jbnd = 1, nbndsub
                   WRITE(iunit, '(2I6,6F10.5,2I4,3ES16.8)') &
                        iki, ikf, xki(:,iki), xkf(:,ikf), ibnd, jbnd, &
                        ABS(edmatf_b(ibnd, jbnd))**2, &
                        REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       IF (ionode .AND. (MOD(iki, MAX(1,nki/10)) == 0 .OR. iki == nki)) THEN
          WRITE(stdout, '(5X,A,I6,A,I6)') '  k_i ', iki, ' / ', nki
          FLUSH(stdout)
       ENDIF
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

    DEALLOCATE(xki, xkf)
  END SUBROUTINE ed_interp_from_file


END MODULE ed_coarse
