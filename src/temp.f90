
SUBROUTINE cal_g_and_det_FT(Q_l, Q_r, D_l, D_r, R_l, R_r, cal_g)
    implicit none
    ! =========================================================================
    ! State-of-the-art DQMC Scale Separation (Stratification) Algorithm
    ! Splits D into D_gt (>1) and D_lt (<=1) to construct an unconditionally
    ! bounded inner matrix. Zero risk of overflow/underflow at extreme \beta.
    ! =========================================================================
    logical, intent(in) :: cal_g
    complex(dp), intent(in) :: Q_l(Ns, Ns), Q_r(Ns, Ns)
    complex(dp), intent(in) :: D_l(Ns), D_r(Ns)
    complex(dp), intent(in) :: R_l(Ns, Ns), R_r(Ns, Ns)
    
    integer :: i, j
    
    ! --- Scale Arrays ---
    complex(dp) :: Dl_gt(Ns), Dl_lt(Ns)
    complex(dp) :: Dr_gt(Ns), Dr_lt(Ns)
    
    ! --- Matrix Buffers ---
    complex(dp) :: Q_r_Q(Ns, Ns)      ! Stores (Q_l * Q_r)^(-1)
    complex(dp) :: Term1(Ns, Ns)      ! Stores (Dr_gt)^(-1) * Q_r_Q * (Dl_gt)^(-1)
    complex(dp) :: Term2(Ns, Ns)      ! Stores Dr_lt * (R_r * R_l) * Dl_lt
    complex(dp) :: M_inner(Ns, Ns)    ! The bounded core matrix
    
    complex(dp) :: R_temp(Ns, Ns), D_temp(Ns)
    complex(dp) :: g_buffer(Ns, Ns), Q_temp(Ns, Ns), RQ_l(Ns, Ns)

    ! =========================================================================
    ! STAGE 1: Split the Scales
    ! =========================================================================
    do i = 1, Ns
      ! Split D_l
      if (abs(D_l(i)) > 1.0d0) then
        Dl_gt(i) = D_l(i)
        Dl_lt(i) = (1.0d0, 0.0d0)
      else
        Dl_gt(i) = (1.0d0, 0.0d0)
        Dl_lt(i) = D_l(i)
      end if
      
      ! Split D_r
      if (abs(D_r(i)) > 1.0d0) then
        Dr_gt(i) = D_r(i)
        Dr_lt(i) = (1.0d0, 0.0d0)
      else
        Dr_gt(i) = (1.0d0, 0.0d0)
        Dr_lt(i) = D_r(i)
      end if
    end do
    print*,'Scale separation complete:'
    print*,'Dl_gt:', Dl_gt
    print*,'Dl_lt:', Dl_lt
    print*,'Dr_gt:', Dr_gt
    print*,'Dr_lt:', Dr_lt
    ! =========================================================================
    ! STAGE 2: Construct the Bounded Inner Matrix M_inner
    ! M_inner = (Dr_gt)^(-1) * (Q_l * Q_r)^(-1) * (Dl_gt)^(-1) + Dr_lt * (R_r * R_l) * Dl_lt
    ! =========================================================================
    
    ! 1. Calculate Q_l^(-1) -> g_buffer
    g_buffer = Q_l
    call inverse(Ns, g_buffer)
    
    ! 2. Calculate Q_r^(-1) -> RQ_l
    RQ_l = Q_r
    call inverse(Ns, RQ_l)
    
    ! 3. Q_r_Q = (Q_l * Q_r)^(-1) = Q_r^(-1) * Q_l^(-1)
    call gemm(Q_r_Q, RQ_l, g_buffer, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    
    ! 4. Construct Term1: Safely scale rows and columns by 1 / D_gt
    do j = 1, Ns
      do i = 1, Ns
        Term1(i, j) = Q_r_Q(i, j) / (Dr_gt(i) * Dl_gt(j))
      end do
    end do
    
    ! 5. Calculate R_r * R_l -> Term2
    call gemm(Term2, R_r, R_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    
    ! 6. Construct Term2: Scale by D_lt
    do j = 1, Ns
      do i = 1, Ns
        Term2(i, j) = Term2(i, j) * Dr_lt(i) * Dl_lt(j)
      end do
    end do
    
    ! 7. Combine into the perfectly bounded M_inner
    do j = 1, Ns
      do i = 1, Ns
        M_inner(i, j) = Term1(i, j) + Term2(i, j)
      end do
    end do

    ! =========================================================================
    ! STAGE 3: Decomposition & Weight Extraction
    ! =========================================================================
    
    ! 8. QR Decompose the stable M_inner: M_inner = Q_M * D_M * R_M
    call qdr(Ns, Ns, M_inner, R_temp, D_temp)
    ! Now: M_inner = Q_M, R_temp = R_M, D_temp = D_M

    ! 9. Calculate Log-Weight
    ! W = det(Q_r_Q^(-1)) * det(Dr_gt) * det(Dl_gt) * det(M_inner)
    !   = [ 1 / det(Q_r_Q) ] * det(Dr_gt) * det(Dl_gt) * [ det(Q_M)*det(D_M)*det(R_M) ]
    ln_cw = -log(det(Ns, Q_r_Q)) + sum(log(Dr_gt)) + sum(log(Dl_gt)) &
            + log(det(Ns, M_inner)) + sum(log(D_temp)) + log(det(Ns, R_temp))

    ! =========================================================================
    ! STAGE 4: Assemble Green's Function
    ! G = Q_l^(-1) * (Dl_gt)^(-1) * R_M^(-1) * D_M^(-1) * Q_M^(-1) * (Dr_gt)^(-1) * Q_r^(-1)
    ! =========================================================================
    
    ! (a) Safely invert R_M
    call inverse(Ns, R_temp) 
    
    ! (b) L_part = Q_l^(-1) * (Dl_gt)^(-1) * R_M^(-1)
    ! Scale columns of Q_l^(-1) (in g_buffer) by 1 / Dl_gt
    do j = 1, Ns
      do i = 1, Ns
        g_buffer(i, j) = g_buffer(i, j) / Dl_gt(j)
      end do
    end do
    ! L_part -> Q_temp
    call gemm(Q_temp, g_buffer, R_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    
    ! (c) Apply D_M^(-1) to columns of L_part
    do j = 1, Ns
      do i = 1, Ns
        Q_temp(i, j) = Q_temp(i, j) / D_temp(j)
      end do
    end do
    
    ! (d) R_part = Q_M^(-1) * (Dr_gt)^(-1) * Q_r^(-1)
    ! Q_M^(-1) is the conjugate transpose of Q_M (in M_inner) -> store in Term1
    do j = 1, Ns
      do i = 1, Ns
        Term1(i, j) = conjg(M_inner(j, i))
      end do
    end do
    
    ! Scale rows of Q_r^(-1) (in RQ_l) by 1 / Dr_gt
    do j = 1, Ns
      do i = 1, Ns
        RQ_l(i, j) = RQ_l(i, j) / Dr_gt(i)
      end do
    end do
    
    ! R_part -> Term2
    call gemm(Term2, Term1, RQ_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    
    ! (e) Final G = L_part * R_part
    call gemm(g_buffer, Q_temp, Term2, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

    ! =========================================================================
    ! STAGE 5: Final Adjustments
    ! =========================================================================
    if(cal_g) then
      g_h(:, :) = - g_buffer(:, :)
      do i = 1, Ns
        g_h(i, i) = g_h(i, i) + 1d0
      end do
    end if

    if(TR_weight) then
      ln_cw = ln_cw + conjg(ln_cw)
    end if
    ln_cw = ln_cw * ncopy
    ln_cw = cmplx(real(ln_cw), mod(aimag(ln_cw), 2*pi), kind=dp)

  END SUBROUTINE cal_g_and_det_FT













  SUBROUTINE cal_g_and_det_FT(Q_l, Q_r, D_l, D_r, R_l, R_r, cal_g)
    implicit none
    ! for FTQMC, g = (1 + B(\tau,0) * B(\beta,\tau))^(-1)
    ! = (1 + Q_r*D_r*R_r * R_l*D_l*Q_l)^(-1)
    ! = (Q_l)^(-1) * [(Q_l Q_r)^(-1) + D_r*(R_r*R_l)*D_l]^(-1) * (Q_r)^(-1)
    ! = (R Q_l)^(-1) * (D)^(-1) * (Q_r Q )^(-1)
    ! where QDR = (Q_l Q_r)^(-1) + D_r * (R_r*R_l) * D_l
    ! ln_cw = log(det(1 +  B(\beta,0))) = log(det(1 + B(\tau,0) * B(\beta,\tau))) = log(det(g^(-1)))
    ! = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D)) for FTQMC
    ! use gemm for matrix multiplication for better performance
    integer :: i,j
    logical,intent(in) :: cal_g
    complex(dp), intent(in) :: Q_l(Ns, Ns), Q_r(Ns, Ns), D_l(Ns), D_r(Ns), R_l(Ns, Ns), R_r(Ns, Ns)
    complex(dp) :: Q_l_Q_r_inv(Ns,Ns), QDR_temp(Ns,Ns), RQ_l(Ns,Ns), Q_r_Q(Ns,Ns)
    complex(dp) :: R_temp(Ns,Ns),D_temp(Ns) ! for QDR calculation
    complex(dp) :: g_buffer(Ns, Ns)
    complex(dp) :: Q_temp(Ns, Ns)
    ln_cw = 0d0
    ! initialization
    QDR_temp = 0d0
    Q_l_Q_r_inv = 0d0
    RQ_l = 0d0
    Q_r_Q = 0d0
    g_buffer = 0d0
    ! Q_l_Q_r_inv = (Q_l Q_r)^(-1)
    call gemm(Q_l_Q_r_inv, Q_l, Q_r, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    call inverse(Ns, Q_l_Q_r_inv)
    
    ! QDR = R_r*R_l
    call gemm(QDR_temp, R_r, R_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

    ! QDR = (Q_l Q_r)^(-1) + D_r*(R_r*R_l)*D_l
    do j = 1, Ns
      do i = 1, Ns
        QDR_temp(i, j) = QDR_temp(i, j) * D_r(i) * D_l(j) + Q_l_Q_r_inv(i, j)
      end do
    end do
    ! do qdr decomposition to get Q, D and R for QDR
    call qdr(Ns, Ns, QDR_temp, R_temp, D_temp)

    ! calculate ln_cw = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D))
    call gemm(RQ_l, R_temp, Q_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    call gemm(Q_r_Q, Q_r, QDR_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    ln_cw = ln_cw + log(det(Ns, Q_r_Q)) + log(det(Ns, RQ_l)) + sum(log(D_temp))
    ! calculate g = (R Q_l)^(-1) * (D)^(-1) * (Q_r Q )^(-1)
    call inverse(Ns, RQ_l)
    call inverse(Ns, Q_r_Q)
    if(cal_g) then
      ! 【修复点】在 gemm 之前，先计算 (R Q_l)^(-1) * D^(-1)
      ! 数学意义：将 RQ_l_inv 的第 j 列除以 D_temp(j)
      ! Fortran 内存优化：内层循环 i 是连续访问的
      do j = 1, Ns
        do i = 1, Ns
          RQ_l(i, j) = RQ_l(i, j) / D_temp(j)
        end do
      end do
      
      ! 然后再执行矩阵乘法： [ (R Q_l)^(-1) * D^(-1) ] * (Q_r Q)^(-1)
      call gemm(g_buffer, RQ_l, Q_r_Q, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
      
      ! 空穴格林函数转换
      g_h(:, :) = - g_buffer(:, :)
      do i = 1, Ns
        g_h(i, i) = g_h(i, i) + 1d0
      end do
    end if

    
    ! make sure weight is in the principal branch
    if(TR_weight) then
      ln_cw = ln_cw + conjg(ln_cw)
    end if
    ln_cw = ln_cw*ncopy
    ln_cw = cmplx(real(ln_cw),mod(aimag(ln_cw),2*pi),kind=dp)

end subroutine





 SUBROUTINE cal_g_and_det_FT(Q_l, Q_r, D_l, D_r, R_l, R_r, cal_g)
    implicit none
    ! for FTQMC, g = (1 + B(\tau,0) * B(\beta,\tau))^(-1)
    ! = (1 + Q_r*D_r*R_r * R_l*D_l*Q_l)^(-1)
    ! = (Q_l)^(-1) * [(Q_l Q_r)^(-1) + D_r*(R_r*R_l)*D_l]^(-1) * (Q_r)^(-1)
    ! = (R Q_l)^(-1) * (D)^(-1) * (Q_r Q )^(-1)
    ! where QDR = (Q_l Q_r)^(-1) + D_r * (R_r*R_l) * D_l
    ! ln_cw = log(det(1 +  B(\beta,0))) = log(det(1 + B(\tau,0) * B(\beta,\tau))) = log(det(g^(-1)))
    ! = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D)) for FTQMC
    ! use gemm for matrix multiplication for better performance
    integer :: i,j
    logical,intent(in) :: cal_g
    complex(dp), intent(in) :: Q_l(Ns, Ns), Q_r(Ns, Ns), D_l(Ns), D_r(Ns), R_l(Ns, Ns), R_r(Ns, Ns)
    complex(dp) :: Ql_temp(Ns, Ns), Qr_temp(Ns, Ns)
    complex(dp) :: Q_l_Q_r(Ns,Ns), R_r_R_l(Ns,Ns),QDR_temp(Ns,Ns), R_Ql(Ns,Ns), Qr_Q(Ns,Ns)
    complex(dp) :: R_temp(Ns,Ns),D_temp(Ns),Q_temp(Ns,Ns) ! for QDR calculation
    complex(dp) :: g_buffer(Ns, Ns)
    ! --- Scale Arrays ---
    complex(dp) :: Dl_gt(Ns), Dl_lt(Ns)
    complex(dp) :: Dr_gt(Ns), Dr_lt(Ns)
    ln_cw = 0d0
    ! initialization
    QDR_temp = 0d0
    Q_l_Q_r = 0d0
    R_r_R_l = 0d0
    R_Ql = 0d0
    Qr_Q = 0d0
    g_buffer = 0d0
    Dl_gt = 0d0
    Dl_lt = 0d0
    Dr_gt = 0d0
    Dr_lt = 0d0
    Q_temp = 0d0
    ! =========================================================================
    ! STAGE 1: Split the Scales
    ! =========================================================================
    do i = 1, Ns
      ! Split D_l
      if (abs(D_l(i)) > 1.0d0) then
        Dl_gt(i) = D_l(i)
        Dl_lt(i) = (1.0d0, 0.0d0)
      else
        Dl_gt(i) = (1.0d0, 0.0d0)
        Dl_lt(i) = D_l(i)
      end if

      ! Split D_r
      if (abs(D_r(i)) > 1.0d0) then
        Dr_gt(i) = D_r(i)
        Dr_lt(i) = (1.0d0, 0.0d0)
      else
        Dr_gt(i) = (1.0d0, 0.0d0)
        Dr_lt(i) = D_r(i)
      end if
    end do

    ! Q_l_Q_r = (Q_l Q_r)^(-1), use that Q_l and Q_r are unitary to avoid matrix inverse
    Qr_temp = conjg(transpose(Q_r))
    Q_l_Q_r = 0d0
    Ql_temp = conjg(transpose(Q_l))
    call gemm(Q_l_Q_r, Qr_temp, Ql_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

    ! QDR = R_r*R_l
    call gemm(R_r_R_l, R_r, R_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

    ! Scale R_r*R_l by Dr_lt and Dl_lt, and scale (Q_l Q_r)^(-1) by Dl_gt and Dr_gt
    do j = 1, Ns
      do i = 1, Ns
        R_r_R_l(i, j) = R_r_R_l(i, j) * Dr_lt(i) * Dl_lt(j)
        Q_l_Q_r(i, j) = Q_l_Q_r(i, j) / Dr_gt(i) / Dl_gt(j)
      end do
    end do

    ! QDR = (Dr_gt)^(-1) * (Q_l Q_r)^(-1) * (Dl_gt)^(-1) + Dr_lt*(R_r*R_l)*Dl_lt
    QDR_temp = R_r_R_l + Q_l_Q_r

    ! QDR_temp = Dr_gt * QDR_temp * Dl_gt
    do j = 1, Ns
      do i = 1, Ns
        QDR_temp(i, j) = QDR_temp(i, j) * Dr_gt(i) * Dl_gt(j)
      end do
    end do
    ! do qdr decomposition to get Q, D and R for QDR
    call qdr(Ns, Ns, QDR_temp, R_temp, D_temp)

    ! calculate ln_cw = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D))
    Qr_temp = Q_r
    Ql_temp = Q_l
    ln_cw = ln_cw + log(det(Ns, Ql_temp)) + log(det(Ns, Qr_temp)) + sum(log(D_temp))  + log(det(Ns,QDR_temp)) + log(det(Ns, R_temp))

    if(cal_g) then
      ! cal g = Ql^(-1) * Dl_gt^(-1) * (R)^(-1) * (D)^(-1) * (Q)^(-1) * Dr_gt^(-1) * Qr^(-1)
      ! Q_r and Q_l are unitary, so their inverse is just their conjugate transpose, no need to do matrix inverse

      call inverse(Ns, R_temp)
      do j = 1, Ns
        do i = 1, Ns
          Ql_temp(i, j) = conjg(Q_l(j, i))
          Qr_temp(i, j) = conjg(Q_r(j, i))
          Q_temp(i, j) = conjg(QDR_temp(j, i)) ! 将 Q_M 的逆存入 Q_temp 以防覆盖
        end do
      end do

      ! g = Q_l^(-1) * R_Temp^(-1) * D_Temp^(-1) * QDR_temp^(-1) * Q_r^(-1)

      call gemm(R_Ql,  Ql_temp,R_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0)) ! R^(-1) * Q_l^(-1)

      call gemm(Qr_Q,  Q_temp, Qr_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0)) ! QDR^(-1) * Q_r^(-1)

      ! cal g_buffer = R_Ql * D_Temp^(-1) * Qr_Q
      do j = 1, Ns
        R_Ql(:, j) = R_Ql(:, j) / D_temp(j)
      end do
      call gemm(g_buffer, R_Ql, Qr_Q, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

      ! 空穴格林函数转换
      g_h(:, :) = - g_buffer(:, :)
      do i = 1, Ns
        g_h(i, i) = g_h(i, i) + 1d0
      end do
    end if

    ! make sure weight is in the principal branch
    if(TR_weight) then
      ln_cw = ln_cw + conjg(ln_cw)
    end if
    ln_cw = ln_cw*ncopy
    ln_cw = cmplx(real(ln_cw),mod(aimag(ln_cw),2*pi),kind=dp)

  end subroutine