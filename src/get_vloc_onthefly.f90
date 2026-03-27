subroutine get_vloc_onthefly(prefix_in, outdir_in, vf_out, is_noncolin, &
                              bxc_out)
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : prefix, tmp_dir
  USE scf,              ONLY : v, vltot
  USE fft_base,         ONLY : dfftp
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,        ONLY : at, alat, celldm, ibrav, omega
  USE edic_mod,         ONLY : V_file

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: prefix_in, outdir_in
  TYPE(V_file), INTENT(INOUT)   :: vf_out
  LOGICAL, INTENT(IN)           :: is_noncolin
  REAL(DP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: bxc_out(:,:)

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  LOGICAL :: wfc_is_collected

  prefix = TRIM(prefix_in)
  tmp_dir = trimcheck(outdir_in)

  wfc_is_collected = .false.
  CALL read_file_new(wfc_is_collected)

  vf_out%nr1 = dfftp%nr1
  vf_out%nr2 = dfftp%nr2
  vf_out%nr3 = dfftp%nr3
  vf_out%nr1x = dfftp%nr1x
  vf_out%nr2x = dfftp%nr2x
  vf_out%nr3x = dfftp%nr3x

  ALLOCATE(vf_out%pot(dfftp%nnr))
  vf_out%pot(:) = v%of_r(:,1) + vltot(:)

  vf_out%nat = nat
  vf_out%ntyp = ntyp
  ALLOCATE(vf_out%ityp(nat))
  ALLOCATE(vf_out%tau(3, nat))
  vf_out%ityp(:) = ityp(1:nat)
  vf_out%tau(:,:) = tau(:, 1:nat)

  vf_out%alat = alat
  vf_out%ibrav = ibrav
  vf_out%celldm = celldm
  vf_out%at = at
  vf_out%omega = omega

  IF (is_noncolin .AND. PRESENT(bxc_out)) THEN
    ALLOCATE(bxc_out(dfftp%nnr, 3))
    bxc_out(:,1) = v%of_r(:,2)
    bxc_out(:,2) = v%of_r(:,3)
    bxc_out(:,3) = v%of_r(:,4)
  ENDIF

end subroutine get_vloc_onthefly
