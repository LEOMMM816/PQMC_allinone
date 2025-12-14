MODULE nrtype
  ! 引入标准环境模块，获取明确的类型定义
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
       int8, int16, int32, &  ! 对应 1, 2, 4 字节整数
       real32, real64         ! 对应 单精度, 双精度
  IMPLICIT NONE
  PRIVATE :: int8, int16, int32, real32, real64
  ! ==========================================
  ! 重构部分：映射到标准类型
  ! ==========================================

  ! Integers:
  ! SELECTED_INT_KIND(9)  -> range > 10^9 -> 32-bit (4 bytes)
  INTEGER, PARAMETER :: I4B = int32
  ! SELECTED_INT_KIND(4)  -> range > 10^4 -> 16-bit (2 bytes)
  INTEGER, PARAMETER :: I2B = int16
  ! SELECTED_INT_KIND(2)  -> range > 10^2 -> 8-bit  (1 byte)
  INTEGER, PARAMETER :: I1B = int8

  ! Reals:
  ! KIND(1.0)   -> Default Real   -> 32-bit (Single)
  INTEGER, PARAMETER :: SP = real32
  ! KIND(1.0D0) -> Default Double -> 64-bit (Double)
  INTEGER, PARAMETER :: DP = real64

  ! Complex:
  ! Fortran中，Complex的Kind值等于其对应的Real的Kind值
  ! 所以 Complex(Kind=real32) 就是单精度复数
  INTEGER, PARAMETER :: SPC = real32
  INTEGER, PARAMETER :: DPC = real64
  
  ! Logical:
  ! ISO_FORTRAN_ENV 标准中通常不直接定义 logical32 等
  ! 保持原样最安全，代表编译器的默认逻辑型
  INTEGER, PARAMETER :: LGT = KIND(.true.)

  ! ==========================================
  ! 以下部分保持原样 (数学常量)
  ! ==========================================
  ! 注意：因为上面定义了 SP 和 DP，这里的 _sp 和 _dp 后缀依然有效
  
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

  ! ==========================================
  ! 以下部分保持原样 (派生类型)
  ! ==========================================

  TYPE sprs2_sp
    INTEGER(I4B) :: n,len
    REAL(SP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_sp

  TYPE sprs2_dp
    INTEGER(I4B) :: n,len
    REAL(DP), DIMENSION(:), POINTER :: val
    INTEGER(I4B), DIMENSION(:), POINTER :: irow
    INTEGER(I4B), DIMENSION(:), POINTER :: jcol
  END TYPE sprs2_dp

END MODULE nrtype
 
  