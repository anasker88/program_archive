*********************************************************************
* �t�@�C����  �FMAT2D.FOR                                           *
* �^�C�g��    �F���`�V�X�e����@�v���O����                          *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� ���p���w��                      *
* �����      �F2003.12.25                                          *
* ����        �FFORTRAN                                             *
*********************************************************************
*�m���e�n                                                           *
*    �{�v���O�����͖{����15.6�߂ɂ���2�����̗��ō쐬���ꂽ���`    *
*    �V�X�e����@��1��ł���D                                      *
*    ��@�́C���ږ@�C�����@�C�����ăN�����t������Ԗ@��p���Ă���D *
*    ��@�A���S���Y���Ƃ��̃T�u���[�`�����͈ȉ��Ɏ����ʂ�ł���D   *
*    -------------------------------------------------------------  *
*    ���ږ@                              �T�u���[�`����             *
*    -------------------------------------------------------------  *
*    LU����@(���I��t)               ---> DLUP                     *
*    Gauss�̏����@(���I��t)          ---> DFGP                     *
*    Gauss�̏����@(���I��)          ---> DFG                      *
*    Gauss-Jordan�@(���I��t)         ---> DGJP                     *
*    Gauss�̏����@(2����(5�_)����     ---> DGBND                    *
*    �o���h�}�g���b�N�X�p: ���I��)                                *
*    -------------------------------------------------------------  *
*    �����@                              �T�u���[�`����             *
*    -------------------------------------------------------------  *
*    JACOBI�@                         ---> HJACOB                   *
*    point-SOR�@                      ---> HSOR                     *
*    point-SOR�@(�o���h�}�g���b�N�X�p)---> SORBND                   *
*    line-SOR�@(�o���h�}�g���b�N�X�p) ---> LBLBND                   *
*    ����: SOR�@�ŁC��=1�Ƃ����Gauss-Seidel�@�ƂȂ�D              *
*    -------------------------------------------------------------  *
*    �N�����t������Ԗ@                  �T�u���[�`����             *
*    -------------------------------------------------------------  *
*    �����c���@                       ---> KCR                      *
*    �����c���@(�o���h�}�g���b�N�X�p) ---> KCRBND                   *
*    Bi-CGSTAB                        ---> KBICG                    *
*    Bi-CGSTAB(�o���h�}�g���b�N�X�p)  ---> KBIBND                   *
*�m�@��@��@�n                                                     *
*    2�����Δ��������� ��^{2}f + ��_{x}f + ��_{y}f = b �ɂ��āC2��*
*  ���x���S�����ߎ���p���ė��U������D�����Ĉȉ��Ɏ����v�Z�̈���l *
*  ���C�ʂ��ԍ�������ꂽ�v�Z�i�q�̔ԍ����̂��̂����ƂȂ�悤��   *
*  ���`�V�X�e��(AX=B)�����C�e��̉�@��p���ĉ����D               *
*        +------------------------+ �i���j                          *
* (NY)5  |  5 | 10 | 15 | 20 | 25 |  f �̒�`�_�͊e�i�q�̒��S�ł��� *
*        +------------------------+  �Ƃ���Df �̒�`�_�ԋ�������� *
*     4  |  4 |  9 | 14 | 19 | 24 |  �e�v�Z�i�q�̕��͂��������=1�� *
*        +------------------------+  ����D                         *
*     3  |  3 |  8 | 13 | 18 | 23 |  �܂��C���E�����Ɋւ��ẮC�v�Z *
*        +------------------------+  �̈�O�̒�`�_�̒l���[���Ƃ��� *
*     2  |  2 |  7 | 12 | 17 | 22 |  �����D                         *
*        +------------------------+  �ʂ��ԍ�:k=(i-1)*NY+j          *
*     1  |  1 |  6 | 11 | 16 | 21 |                                 *
*        +------------------------+                                 *
*      j   i=1     2    3    4    5(NX)                             *
*  x,y�������ꂼ��5�������C�ʂ��ԍ���k=1����25�܂ł��Cf{i,j}=k��  *
*  ���ƂȂ�悤�ɂ���D                                             *
* [�@���@���@��@]                                                  *
* ( 1 )�@���@�U�@��                                                 *
*   ��^{2}f + ��_{x}f + ��_{y}f = b                                 *
*   ===> ( f_{i-1,j  }-2f_{i,j}+f_{i+1,j  } ) / ��^{2}              *
*       +( f_{i  ,j-1}-2f_{i,j}+f_{i  ,j+1} ) / ��^{2}              *
*       +( f_{i+1,j  }-f_{i-1,j  } ) / (2��)                        *
*       +( f_{i  ,j+1}-f_{i  ,j-1} ) / (2��)                        *
*       = b_{k}   where k=(i-1)*NY + j                              *
*   ===> (  1/��^{2} - 1/(2��)  ) f_{i-1,j  }                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~-> AT1(k) = A(i-1,j)               *
*       +( -2/��^{2} - 2/��^{2} ) f_{i  ,j  }                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~-> AT3(k) = A(i,j)                 *
*       +(  1/��^{2} + 1/(2��)  ) f_{i+1,j  }                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~-> AT5(k) = A(i+1,j)               *
*       +(  1/��^{2} - 1/(2��) )  f_{i  ,j-1}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~-> AT2(k) = A(i,j-1)               *
*       +(  1/��^{2} + 1/(2��) )  f_{i  ,j+1}                       *
*        ~~~~~~~~~~~~~~~~~~~~~~~~-> AT4(k) = A(i,j+1)               *
*       = b_{k}   where k=(i-1)*NY + j  : ���ӁF��=1�Ƃ���          *
*   ===> Af = B ===> AX = B                                         *
* ( 2 )�@�}�g���b�N�X�̐ݒ�                                         *
*  2-1 A�̐ݒ�                                                      *
*    A��65�s����73�s�ɂ���悤�ɐݒ肷��DAT1,AT2,AT3,AT4,AT5 ��    *
*    �o���h�}�g���b�N�X�p�̉�@�̂��߂̂��̂ł���D����ȊO�́C�[�� *
*    �̗v�f���܂ޑS(�t��)�}�g���b�N�X�������ꍇ�́CA��p����D      *
*  2-2 B�̐ݒ�                                                      *
*    65�s����73�s�ɂ����āC���ƂȂ�f�̒l�����������̂�B�Ƃ���΂� *
*    ���D��̓I�ɂ́C�ȉ��̂悤�ɂ��ď��Ƀ}�g���b�N�X�����߂Ă䂭�D *
*    k= 1 (i=1,j=1)->f_{i-1,j  }:�̈�O�ƂȂ�̂�0�Ƃ���            *
*                    f_{i  ,j  }:(i-1)*NY+j����                   *
*                    f_{i+1,j  }:((i+1)-1)*NY+j����               *
*                    f_{i  ,j-1}:�̈�O�ƂȂ�̂�0�Ƃ���            *
*                    f_{i  ,j+1}:(i-1)*NY+(j+1)����               *
*                    b_{k}:�ȏ��f�̒l�ƃ�(=1)���狁�܂�            *
*     k= 2 (i=1,j=2), k= 3 (i=1,j=3), k= 4 (i=1,j=4),               *
*     k= 5 (i=1,j=5), k= 6 (i=2,j=1), k= 7 (i=2,j=2),               *
*     ...... ,        k=25 (i=5,j=5)                                *
*   i=3, j=3 (k=13�ɂ�����CAT1,AT2,AT3,AT4,AT5��A�̊֌W            *
*        +------------------------+                                 *
* (NY)5  |    |    |    |    |    |                                 *
*        +------------------------+                                 *
*     4  |    |    |AT4 |    |    |    k=(i-1)*NY+j                 *
*        +------------------------+                                 *
*     3  |    | AT1|AT3 |AT5 |    |    AT1(k) -> A(k,k-NY)          *
*        +---------+--------------+    AT2(k) -> A(k,k- 1)          *
*     2  |    |    |AT2 |    |    |    AT3(k) -> A(k,k   )          *
*        +------------------------+    AT4(k) -> A(k,k+ 1)          *
*     1  |    |    |    |    |    |    AT5(k) -> A(k,k+NY)          *
*        +------------------------+    B(k)   -> B(k)               *
*     j   i=1     2    3    4    5(NX)                              *
*********************************************************************
      PROGRAM MAT2D
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
* �p�����[�^�ϐ�(NX0��0�����Ă���̂̓T�u���[�`�������ɂ�PARAMETER
* �ϐ��ł���NX0�𒼐ړn���Ȃ����߁D���NX=NX0�Ƒ���������Ă����n��
* �Ă���D)  NX:x�����i�q��, NY:y�����i�q��, NE:�S�i�q��=NX*NY
      PARAMETER ( NX0=5, NY0=5, NE0=25 )
      DIMENSION A(NE0,NE0),B(NE0),X1(NE0),X2(NX0,NY0)
* �o���h�}�g���b�N�X�p�z��̒�`
      DIMENSION AT1(NE0),AT2(NE0),AT3(NE0),AT4(NE0),AT5(NE0)
* ��Ɨp�z��
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),W6(NE0),
     $          W7(NE0),W8(NE0),AGB(-NY0:NY0,NE0)
      DIMENSION IPV(NE0)
* ��L���U���������̃�=1��ݒ�
      DELTA = 1.0D0
* PARAMETER�ϐ��̒l��SUBROUTINE�����ɓn�����߂ɑ��
      NX = NX0
      NY = NY0
      NE = NE0
* �߂�_ : ��@�̕ύX���̖߂�_
  700 CONTINUE
* A,B ��ݒ肷�邽�߂ɂ��炩����X2�̒l���i�q�ԍ��Ɠ����ɂȂ�悤�ɂ���
      I = 1
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          X2(IX,IY) = DBLE(I)
          I =I + 1
   20   CONTINUE
   10 CONTINUE
* X2 �ȊO�̔z��̃[���N���A
      DO 30 I = 1,NE
        AT1(I) = 0.0D0
        AT2(I) = 0.0D0
        AT3(I) = 0.0D0
        AT4(I) = 0.0D0
        AT5(I) = 0.0D0
        B(I) = 0.0D0
        X1(I) = 0.0D0
        W1(I)=0.0D0
        W2(I)=0.0D0
        W3(I)=0.0D0
        W4(I)=0.0D0
        W5(I)=0.0D0
        W6(I)=0.0D0
        W7(I)=0.0D0
        W8(I)=0.0D0
        DO 40 II = 1,NE
          A(I,II) = 0.0D0
   40   CONTINUE
        DO 45 III = -NY,NY
          AGB(III,I) = 0.0D0
   45   CONTINUE
   30 CONTINUE
* X_{i,j}=k=(i-1)*NY+j �ƂȂ�悤�� A,B �����߂�
      I = 1
      DO 50 IX = 1,NX
        DO 60 IY = 1,NY
*         f(i-1,j)���v�Z�̈�O�i�����j�̋��E�����̏���
          IF (IX.EQ.1) THEN
*           f(0,j)=0�Ƃ��ď�������
            AT1(I) = 0.0D0
*           AT1*f(0,j)��B1�Ƃ���
            B1 = 0.0D0
*         f(i-1,j)���v�Z�̈���̏ꍇ
          ELSE
            AT1(I) = 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA)
            A(I,I-NY) = AT1(I)
            B1 = X2(IX-1,IY  )*( 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA) )
          END IF
*         f(i,j)�͏�Ɍv�Z�̈��
          AT3(I) = -2.0D0/DELTA**2 - 2.0D0/DELTA**2
          A(I,I) = AT3(I)
*         AT3*f(i,j)��B3�Ƃ���
          B3 = X2(IX  ,IY  )*( -2.0D0/DELTA**2 - 2.0D0/DELTA**2 )
*         f(i+1,j)���v�Z�̈�O�i�E���j�̋��E�����̏���
          IF (IX.EQ.NX) THEN
*           f(NX,j)=0�Ƃ��ď���
            AT5(I) = 0.0D0
*           AT5*f(NX,j)=B5�Ƃ���
            B5 = 0.0D0
*         f(i+1,j)���v�Z�̈���̏ꍇ
          ELSE
            AT5(I) = 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA)
            A(I,I+NY) = AT5(I)
            B5 = X2(IX+1,IY  )*( 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA) )
          END IF
*         f(i,j-1)���v�Z�̈�O�i�����j�̋��E�����̏���
          IF (IY.EQ.1) THEN
*           f(i,j-1)=0�Ƃ��ď���
            AT2(I) = 0.0D0
*           AT2*f(i,j-1)=B2�Ƃ���
            B2 = 0.0D0
*         f(i,j-1)���v�Z�̈���̏ꍇ
          ELSE
            AT2(I) = 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA)
            A(I,I-1) = AT2(I)
            B2 = X2(IX  ,IY-1)*( 1.0D0/DELTA**2 - 1.0D0/(2.0D0*DELTA) )
          END IF
*         f(i,j+1)���v�Z�̈�O�i�㑤�j�̋��E�����̏���
          IF (IY.EQ.NY) THEN
*           f(i,NY)=0�Ƃ��ď���
            AT4(I) = 0.0D0
*           AT4*f(i,NY)=B4�Ƃ���
            B4 = 0.0D0
*         f(i,j+1)���v�Z�̈���̏ꍇ
          ELSE
            AT4(I) = 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA)
            A(I,I+1) = AT4(I)
            B4 = X2(IX  ,IY+1)*( 1.0D0/DELTA**2 + 1.0D0/(2.0D0*DELTA) )
          END IF
*         B�̌v�Z
          B(I) = B1 + B2 + B3 + B4 + B5
          I = I + 1
   60   CONTINUE
   50 CONTINUE
* �߂�_ : ��@�̑I���������ȂƂ��̖߂�_
  710 CONTINUE
* �e�T�u���[�`����p���Čv�Z����
      WRITE (6,2000)
 2000 FORMAT ( '  1 : LU Decomposition '/,
     $         '  2 : Gauss Elimination '/,
     $         '  3 : Gauss Elimination without Pivotting '/,
     $         '  4 : Gauss Jordan '/,
     $         '  5 : Gauss Elimination for Band Matrix '/,
     $         '  6 : Jacobi '/,
     $         '  7 : point-SOR '/,
     $         '  8 : point-SOR for Band matrix '/,
     $         '  9 : line-SOR for Band Matrix '/,
     $         ' 10 : Conjugate Residual '/,
     $         ' 11 : Conjugate Residual for Band Matrix '/,
     $         ' 12 : Bi-CGSTAB '/,
     $         ' 13 : Bi-CGSTAB for Band Matrix '/,
     $         ' 14 : END '/,
     $         ' Input Number ---> ')
      READ (5,*) IJ
*   ���ږ@
      IF (IJ.EQ.1) THEN
        CALL DLU (A,B,X1,NE,
     $            W1,IPV)
      ELSE IF (IJ.EQ.2) THEN
        CALL DFGP (A,B,X1,NE,
     $             IPV)
      ELSE IF (IJ.EQ.3) THEN
        CALL DFG (A,B,X1,NE)
      ELSE IF (IJ.EQ.4) THEN
        CALL DGJP (A,B,X1,NE,
     $             IPV)
      ELSE IF (IJ.EQ.5) THEN
        CALL DGBND (AT1,AT2,AT3,AT4,AT5,B,X1,NY,NE,
     $              AGB)
*   �����@
      ELSE IF (IJ.EQ.6) THEN
        CALL HJACOB (A,B,X1,NE,
     $               W1)
      ELSE IF (IJ.EQ.7) THEN
        CALL HSOR (A,B,X1,NE)
      ELSE IF (IJ.EQ.8) THEN
        CALL SORBND (AT1,AT2,AT3,AT4,AT5,B,X1,NY,NE)
      ELSE IF (IJ.EQ.9) THEN
        CALL LBLBND (AT1,AT2,AT3,AT4,AT5,B,X1,NX,NY,NE,
     $               W1,W2,W3,W4,W5,W6,W7,W8)
*   �N�����t������Ԗ@
      ELSE IF (IJ.EQ.10) THEN
        CALL KCR (A,B,X1,NE,
     $            W1,W2,W3,W4)
      ELSE IF (IJ.EQ.11) THEN
        CALL KCRBND (AT1,AT2,AT3,AT4,AT5,B,X1,NY,NE,
     $               W1,W2,W3,W4)
      ELSE IF (IJ.EQ.12) THEN
        CALL KBICG (A,B,X1,NE,
     $              W1,W2,W3,W4,W5,W6)
      ELSE IF (IJ.EQ.13) THEN
        CALL KBIBND (AT1,AT2,AT3,AT4,AT5,B,X1,NY,NE,
     $               W1,W2,W3,W4,W5,W6)
      ELSE IF (IJ.EQ.14) THEN
        GO TO 900
      ELSE
        GO TO 710
      END IF
      GO TO 700
  900 CONTINUE
* �v�Z�I��
      STOP
      END
*********************************************************************
*                  ���ږ@ : ��@ 1 - 5
*********************************************************************
*    LU�����ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��    *
*    (�������I��t)        AX=B                                     *
*    LU�����ł͎��̂悤��2�i�K�ŉ������߂�                          *
*    ��1�i�K : LY = B ���� Y �����߂�                               *
*    ��2�i�K : UX = Y ���� X �����߂�                               *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XX                      *
*            ����͑�2�i�K�� UX=Y �̉� X                            *
*            (���̃T�u���[�`���ł�������߂�)                       *
*    Y(NE) : LY=B �̉� Y                                            *
*    IPV(NE) : �������I��p�̔z��                                   *
*             �S���� NE ��̏���������s�����C�e�����i�K�ōs�����Z  *
*             �ɂ��ĕ��ꂪ�ł��傫���Ȃ�悤�ȗ�̔ԍ�������      *
*********************************************************************
      SUBROUTINE DLU (A,B,X,NE,
     $                Y,IPV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE),Y(NE)
      DIMENSION IPV(NE)
* �z��̃[���N���A
      DO 10 I = 1,NE
        Y(I)  = 0.0D0
   10 CONTINUE
* �������I��p�z��̏����ݒ�
      DO 20 I = 1, NE
        IPV(I) = I
   20 CONTINUE
* �������I���ɂ����ِ��̔���l
      EPS = 1.0D-50
      DO 30 K = 1, NE
        L = K
        APV = ABS( A(IPV(L),K) )
        DO 40 I = K+1, NE
*         �������I��
          IF ( ABS( A(IPV(I),K) ). GT. APV ) THEN
            L = I
            APV = ABS( A(IPV(L),K) )
          END IF
   40   CONTINUE
*       �������I�����s���������悢�Ɣ��f���������ւ����s��
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       �������I�����s���Ă��v�Z�s�\�ȂƂ� : �s��͓���(singular)
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K = ', I4)
          STOP
        END IF
* U �����߂�
        A(IPV(K),K) = 1.0D0 / A(IPV(K),K)
        DO 50 I = K+1, NE
          A(IPV(I),K) = A(IPV(I),K) * A(IPV(K),K)
          DO 60 J = K+1, NE
            A(IPV(I),J) = A(IPV(I),J) - A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
   50   CONTINUE
   30 CONTINUE
* ��1�i�K LY = B ������
      Y(1) = B(IPV(1))
      DO 70 I = 2, NE
        T = B(IPV(I))
        DO 80 J = 1, I-1
          T = T - A(IPV(I),J) * Y(J)
   80   CONTINUE
        Y(I) = T
   70 CONTINUE
* ��2�i�K UX = Y ������
      X(NE) = Y(NE) * A(IPV(NE),NE)
      DO 90 I = NE-1, 1, -1
        T = Y(I)
        DO 100 J = I+1, NE
          T = T - A(IPV(I),J) * X(J)
  100   CONTINUE
        X(I) = T * A(IPV(I),I)
   90 CONTINUE
* ���ʏo��
      DO 110 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_1_( ',I3,' ) = ',1PD13.5)
  110 CONTINUE
      RETURN
      END
*********************************************************************
*    Gauss�̏����@(�������I��t)�ɂ���Ώ̍s��A���܂ސ��`�V�X�e�� *
*    ��@�T�u���[�`��    AX=B                                       *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XX                      *
*            (���̃T�u���[�`���ł�������߂�)                       *
*    IPV(NE) : �������I��p�̔z��                                   *
*             �S���� NE ��̏���������s�����C�e�����i�K�ōs�����Z  *
*             �ɂ��ĕ��ꂪ�ł��傫���Ȃ�悤�ȗ�̔ԍ�������      *
*********************************************************************
      SUBROUTINE DFGP (A,B,X,NE,
     $                 IPV)
      IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
      DIMENSION IPV(NE)
* �������I���ɂ����ِ��̔���l
      EPS = 1.0D-50
* �������I��p�z��̏����ݒ�
      DO 10 I = 1, NE
        IPV(I) = I
   10 CONTINUE
      DO 20 K = 1, NE
        L = K
        APV = ABS( A(IPV(L),K) )
*       �������I��
        DO 30 I = K+1, NE
          IF ( ABS( A(IPV(I),K) ). GT. APV ) THEN
            L = I
            APV = ABS( A(IPV(L),K) )
          END IF
   30   CONTINUE
*       �������I�����s���������悢�Ɣ��f���������ւ����s��
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       �������I�����s���Ă��v�Z�s�\�ȂƂ� : �s��͓���
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K= ', I4)
          STOP
        END IF
* �O�i����
        DO 40 I = K+1, NE
          A(IPV(K),I)=A(IPV(K),I)/A(IPV(K),K)
   40   CONTINUE
        B(IPV(K)) = B(IPV(K))/A(IPV(K),K)
        DO 50 I = K+1, NE
          DO 60 J = K+1, NE
            A(IPV(I),J)=A(IPV(I),J)-A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
          B(IPV(I))=B(IPV(I))-A(IPV(I),K)*B(IPV(K)) 
   50   CONTINUE
   20 CONTINUE
* ��ޑ��
      X(NE) = B(IPV(NE))
      DO 70 I = NE-1, 1, -1
        T = B(IPV(I))
        DO 80 J = I+1, NE
          T = T - A( IPV(I), J ) * X(J)
   80   CONTINUE
        X(I) = T
   70 CONTINUE
* ���ʏo��
      DO 90 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_2_( ',I3,' ) = ',1PD13.5)
   90 CONTINUE
      RETURN
      END
*********************************************************************
*    Gauss�̏����@(�������I��)�ɂ���Ώ̍s��A���܂ސ��`�V�X�e�� *
*    ��@�T�u���[�`��    AX=B                                       *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XX ---> ��������߂�    *
*********************************************************************
      SUBROUTINE DFG (A,B,X,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* ���ِ��̔���l
      EPS = 1.0D-50
* �O�i����
      DO 10 K = 1, NE
*       ���Z���s���W�����������v�Z�s�\�ȂƂ� : �s��͓���
        IF ( ABS( A(K, K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K = ', I4)
          STOP
        END IF
        DO 20 I = K+1, NE
          A(K,I)=A(K,I)/A(K,K)
   20   CONTINUE
        B(K) = B(K)/A(K,K)
        DO 30 I = K+1, NE
          DO 40 J = K+1, NE
            A(I,J)=A(I,J)-A(I,K)*A(K,J)
   40     CONTINUE
          B(I)=B(I)-A(I,K)*B(K) 
   30   CONTINUE
   10 CONTINUE
* ��ޑ��
      X(NE) = B(NE)
      DO 50 I = NE-1, 1, -1
        T = B(I)
        DO 60 J = I+1, NE
          T = T - A( I, J ) * X(J)
   60   CONTINUE
        X(I) = T
   50 CONTINUE
* ���ʏo��
      DO 70 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_3_( ',I3,' ) = ',1PD13.5)
   70 CONTINUE
      RETURN
      END
*********************************************************************
*    Gauss-Jordan�@(�������I��t)�ɂ���Ώ̍s��A���܂ސ��`�V�X�e��*
*    ��@�T�u���[�`��    AX=B                                       *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XX ---> ��������߂�    *
*    IPV(NE) : �������I��p�̔z��                                   *
*             �S���� NE ��̏���������s�����C�e�����i�K�ōs�����Z  *
*             �ɂ��ĕ��ꂪ�ł��傫���Ȃ�悤�ȗ�̔ԍ�������      *
*********************************************************************
      SUBROUTINE DGJP (A,B,X,NE,
     $                 IPV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
      DIMENSION IPV(NE)
* �������I���ɂ����ِ��̔���l
      EPS = 1.0D-50
* �������I��p�z��̏����ݒ�
      DO 10 I = 1, NE
        IPV(I) = I
   10 CONTINUE
      DO 20 K = 1, NE
        L = K
        AIPV = ABS( A(IPV(L),K) )
*       �������I��
        DO 30 I = K+1, NE
          IF ( ABS( A(IPV(I),K) ). GT. AIPV ) THEN
            L = I
            AIPV = ABS( A(IPV(L),K) )
          END IF
   30   CONTINUE
*       �������I�����s���������悢�Ɣ��f���������ւ����s��
        IF (L.NE.K) THEN
          IPVEX = IPV(K)
          IPV(K) = IPV(L)
          IPV(L) = IPVEX
        END IF
*       �������I�����s���Ă��v�Z�s�\�ȂƂ�
        IF ( ABS( A(IPV(K),K) ). LE. EPS ) THEN
          WRITE (*,2000) K
 2000     FORMAT (' Matrix is singular at K= ', I4)
          STOP
        END IF
* �����ߒ�
        DO 40 I = K+1, NE
          A(IPV(K),I)=A(IPV(K),I)/A(IPV(K),K)
   40   CONTINUE
        B(IPV(K)) = B(IPV(K))/A(IPV(K),K)
        DO 50 I = 1, NE
          IF (I.EQ.K) GO TO 50
          DO 60 J = K+1, NE
            A(IPV(I),J)=A(IPV(I),J)-A(IPV(I),K)*A(IPV(K),J)
   60     CONTINUE
          B(IPV(I))=B(IPV(I))-A(IPV(I),K)*B(IPV(K)) 
   50   CONTINUE
   20 CONTINUE
* ��̏����ߒ����I������ƁCIX=B�ƂȂ���͂������܂�
      DO 70 I = 1, NE
        X(I) = B(IPV(I))
   70 CONTINUE
* ���ʏo��
      DO 80 I = 1,NE
        WRITE(6,2010) I,X(I)
 2010   FORMAT(' X_4_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*    2����(5�_)�����p�o���h�}�g���b�N�X��@�T�u���[�`�� AX=B        *
*    2����(5�_)�����ɂē���ꂽ�K���I��Ώ̍s��A���܂񂾐��`�V�X�e��*
*    ��Gauss�̏����@��p���ĉ����T�u���[�`���D(�������I��)        *
*    ������@�̂��߂Ƀo���h�}�g���b�N�X�p�ɂ��Ă���D               *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE)                        *
*             i=3, j=3 (k=13)�ɂ�����CAT1,AT2,AT3,AT4,AT5          *
*                +------------------------+                         *
*         (NY)5  |    |    |    |    |    |                         *
*                +------------------------+                         *
*             4  |    |    |AT4 |    |    |    k=(I-1)*NY+J         *
*                +------------------------+                         *
*             3  |    | AT1|AT3 |AT5 |    |    AT1(k)->AT(-NY,k)    *
*                +------------------------+    AT2(k)->AT(- 1,k)    *
*             2  |    |    |AT2 |    |    |    AT3(k)->AT(  0,k)    *
*                +------------------------+    AT4(k)->AT(+ 1,k)    *
*             1  |    |    |    |    |    |    AT5(k)->AT(+NY,k)    *
*                +------------------------+    B(k)                 *
*     j   i=1     2    3    4    5(NX)                              *
*    AT(-NY:NY,NE) : 2���������ߎ��ɂ��K���I��Ώ̍s��            *
*          -NY:NY->��̐}�ōl�����k=13�̂Ƃ�AT1����AT5�̒ʂ��ԍ��� *
*                  AT1�� k"-NY" : AT2�� k"-1"   : AT3�� k"+0"       *
*                  AT4�� k"+1"  : AT5�� k"+NY"                      *
*                  �ƂȂ�D����"��"�ň͂܂ꂽ�l��-NY:NY�ł���D     *
*                  �Ȃ��C����ȊO��AT((-4,-3,-2,2,3,4),NE)�̓[��    *
*                  �ƂȂ�                                           *
*          NE->��q�̂悤��k�͑S����1����NE                         *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XX ---> ��������߂�    *
*********************************************************************
      SUBROUTINE DGBND (AT1,AT2,AT3,AT4,AT5,B,X,NY,NE,
     $                  AT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE)
      DIMENSION AT(-NY:NY,NE),B(NE),X(NE)
* �}�g���b�N�XAT�̃[���N���A
      DO 10 INE = 1,NE
        DO 20 I = -NY,NY
          AT(I,INE) = 0.0D0
   20   CONTINUE
   10 CONTINUE
* AT1����AT5�� AT �Ɋi�[
      DO 30 INE = 1,NE
        AT(-NY,INE) = AT1(INE)
        AT( -1,INE) = AT2(INE)
        AT(  0,INE) = AT3(INE)
        AT(  1,INE) = AT4(INE)
        AT( NY,INE) = AT5(INE)
   30 CONTINUE
* �O�i����
      DO 40 I = 1,NE-1
        IF ( I.LE.NE-NY ) THEN
          DO 50 J = 1,NY
            AA = AT(-J,I+J)/AT(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 60 K = -J+1,NY-J
              AT(K,I+J) = AT(K,I+J)-AT(N,I)*AA
              N = N + 1
   60       CONTINUE
   50     CONTINUE
        ELSE
          DO 70 J = 1,NE-I
            AA = AT(-J,I+J)/AT(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 80 K = -J+1,NE-I-J
              AT(K,I+J) = AT(K,I+J)-AT(N,I)*AA
              N = N + 1
   80       CONTINUE
   70     CONTINUE
        END IF
   40 CONTINUE
*�W���s��̓��ِ��𔻒�
      IF ( DABS(AT(0,NE)).LE.1.0D-20 ) THEN
        WRITE (6,*) ' Matrix is singular : |A(0,NE)| < 1E-20 '
      END IF
* ��ޑ��
      X(NE) = B(NE) / AT(0,NE)
      DO 90 I = NE-1,1,-1
        S = 0.0D0
        IF ( I.GT.NE-NY ) THEN
          DO 100 N = 1,NE-I
            S = S + AT(N,I)* X(I+N)
  100     CONTINUE
          X(I) = ( B(I)-S ) / AT(0,I)
        ELSE
          DO 110 N = 1,NY
            S = S + AT(N,I) * X(I+N)
  110     CONTINUE
          X(I) = ( B(I)-S ) / AT(0,I)
        END IF
   90 CONTINUE
* ���ʏo��
      DO 120 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_5_( ',I3,' ) = ',1PD13.5)
  120 CONTINUE
      RETURN
      END
*********************************************************************
*                  �����@ : ��@ 6 - 9
*********************************************************************
*    Jacobi�@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��  *
*                        AX=B                                       *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : ���l ; XN(NE) : �V�l ---> ��������߂�                 *
*********************************************************************
      SUBROUTINE HJACOB (A,B,XN,NE,
     $                   X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE),XN(NE)
* �ő�J��Ԃ���(NITR)�Ǝ�������l(EITR)�̐ݒ�
      NITR = 200
      EITR = 1.0D-9
* X �̃[���N���A
      DO 10 I = 1,NE
        X(I) = 0.0D0
   10 CONTINUE
* B��2��m�����̌v�Z
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* ITR : �J��Ԃ��񐔂̃J�E���^
      ITR = 0
* �J��Ԃ��̂��߂̖߂�_
  700 CONTINUE
* ��ɓ���ꂽ�l�����l�ɐݒ肵����
      ITR = ITR + 1
      DO 30 I = 1,NE
        X(I) = XN(I)
   30 CONTINUE
* �V�l�̌v�Z
      DO 40 I = 1, NE
        SUM = 0.0D0
        DO 50 J = 1, NE
          IF (I.EQ.J) GO TO 50
          SUM = SUM + A(I,J)*X(J)
   50   CONTINUE
        XN(I) = ( B(I)-SUM )/A(I,I)
   40 CONTINUE
* ��������̂��߂̃m�����̌v�Z
*     RNORM : �c���x�N�g��(R=AX-B)��2��m����
      RNORM = 0.0D0
      DO 60 I = 1,NE
        SUM = 0.0D0
        DO 70 J = 1,NE
          SUM = SUM + A(I,J)*XN(J)
   70   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   60 CONTINUE
*     �c���̌v�Z ZANSA = || R || / || B ||
      ZANSA = DSQRT(RNORM/BNORM)
* ��������
      IF (ITR.LE.NITR) THEN
        IF (ZANSA.GT.EITR) THEN
          GO TO 700
        ELSE
          WRITE (6,*) 'Converged : Total ITR = ',ITR
        END IF
      ELSE
        WRITE (6,*) ' Not converged ! '
      END IF
* ���ʏo��
      DO 80 I = 1,NE
        WRITE(6,2000) I,XN(I)
 2000   FORMAT(' X_6_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*   point-SOR�@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*    (�ɘa�W��OMG��1�Ƃ����Gauss-Seidel�@)  AX=B                   *
*    NX : x�����i�q������; NY : y�����i�q������; NE : ���i�q�_��    *
*    A(NE,NE) : AX=B �̃}�g���b�N�X A                               *
*    B(NE) : ���`�V�X�e�� AX=B �̃}�g���b�N�XB                      *
*    X(NE) : �}�g���b�N�X X ��1�����z��                             *
*********************************************************************
      SUBROUTINE HSOR (A,B,X,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* �ő�J��Ԃ���(NITR)�Ǝ�������l(EITR)�̐ݒ�
      NITR = 200
      EITR = 1.0D-9
* �ɘa�W���̐ݒ� ( OMG=1.0D0�Ȃ�Gauss-Seidel�@ )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* B��2��m�����̌v�Z
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* ITR : �J��Ԃ��񐔂̃J�E���^
      ITR = 0
* �J��Ԃ��̂��߂̖߂�_
  700 CONTINUE
      ITR = ITR + 1
      DO 20 I = 1, NE
        XOLD = X(I)
        SUM = 0.0D0
        DO 30 J = 1, NE
          IF (I.EQ.J) GO TO 30
          SUM = SUM + A(I,J)*X(J)
   30   CONTINUE
        XNEW = ( B(I)-SUM )/A(I,I)
        X(I) = XOLD + OMG * ( XNEW - XOLD )
   20 CONTINUE
* ��������̂��߂̃m�����̌v�Z
*     RNORM : �c���x�N�g��(R=AX-B)��2��m����
      RNORM = 0.0D0
      DO 40 I = 1,NE
        SUM = 0.0D0
        DO 50 J = 1,NE
          SUM = SUM + A(I,J)*X(J)
   50   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   40 CONTINUE
*     �c���̌v�Z ZANSA = || R || / || B ||
      ZANSA = DSQRT(RNORM/BNORM)
* ��������
      IF (ITR.LE.NITR) THEN
        IF (ZANSA.GT.EITR) THEN
          GO TO 700
        ELSE
          WRITE (6,*) 'Converged : Total ITR = ',ITR
        END IF
      ELSE
        WRITE (6,*) ' Not converged ! '
      END IF
* ���ʏo��
      DO 60 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_7_( ',I3,' ) = ',1PD13.5)
   60 CONTINUE
      RETURN
      END
*********************************************************************
*  point-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*  (�o���h�}�g���b�N�X�p)                                           *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*�m�ϐ��̐����n                                                     *
*    NX:x�����i�q������; NY:y�����i�q������; NE:���i�q�_��= NX*NY   *
*    NITR : �ő唽���� ; EITR : ��������l                        *
*********************************************************************
      SUBROUTINE SORBND (A1,A2,A3,A4,A5,B,X,NY,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),B(NE),X(NE)
* �ő�J��Ԃ���(NITR)�Ǝ�������l(EITR)�̐ݒ�
      NITR = 200
      EITR = 1.0D-9
* �ɘa�W���̐ݒ� ( OMG=1.0D0�Ȃ�Gauss-Seidel�@ )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* B��2��m�����̌v�Z
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* ITR : �J��Ԃ��񐔂̃J�E���^
      DO 20 ITR = 1,NITR
*       RNORM : �c���x�N�g��(R=AX-B)��2��m����
        RNORM = 0.0D0
*       A3,A4,A5�͈̔�
        I=1
          XOLD = X(I)
          SUM =  A4(I)*X(I+1)+A5(I)*X(I+NY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
*       A2,A3,A4,A5�͈̔�
        DO 30 I=2,NY
          XOLD = X(I)
          SUM =  A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   30   CONTINUE
*       A1 - A5 �͈̔�
        DO 40 I=NY+1,NE-NY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   40   CONTINUE
*       A1 - A4 �͈̔�
        DO 50 I=NE-NY+1,NE-1
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
   50   CONTINUE
*       A1 - A3 �͈̔�
        I=NE
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
*       RNORM : �c���x�N�g��(R=AX-B)��2��m����
        RNORM= 0.0D0
        DO 60 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NE-NY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   60   CONTINUE
*       �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ��������
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) 'Converged : Total ITR = ',ITR
          GO TO 700
        END IF
   20 CONTINUE
* NITR�܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not converged ! '
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  700 CONTINUE
* ���ʏo��
      DO 70 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_8_( ',I3,' ) = ',1PD13.5)
   70 CONTINUE
      RETURN
      END
*********************************************************************
*  line-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`�� *
*  (�o���h�}�g���b�N�X�p)                                           *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*    A_{i,j} ---> AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE)           *
*    B -> BX(NE) : ���m�x�N�g��                                     *
*    X -> XN(NE) : ���m�x�N�g�� ---> ��������߂�                   *
* [�g�[�}�X�@�̂��߂̌W���s��]                                      *
*    x,y�����ɉA�I�ɗ��U�����ꂽ���ʂ��ȉ��̂悤�ɕ\���D            *
*  |B(1) C(1) 0    0 ...             |X(1)   |  |D(1)   |           *
*  |A(2) B(2) C(2) 0 ...             |X(2)   |  |D(2)   |           *
*  |0    A(3) B(3) C(3) 0 ...        |X(3)   |= |D(3)   |           *
*  |              ...                |...    |  |...    |           *
*  |0   ...   A(NE-1) B(NE-1) C(NE-1)|X(NE-1)|  |D(NE-1)|           *
*  |0    0    0 ....  A(NE  ) B(NE)  |X(NE)  |  |D(NE)  |           *
*�m�ϐ��̐����n                                                     *
*    NX:x�����i�q������; NY:y�����i�q������; NE:���i�q�_��=NX*NY    *
*    NITR : �ő唽����; EITR : ��������l                         *
*    OMG : �ɘa�W���D1.0�ŏ\���D                                    *
*           ���ӁFPoint-SOR�ƈقȂ�C���܂�傫����������Ɣ��U���� *
* [�z��̐���]                                                      *
* XN...�e�����ւ̑|�����X(�ԍ��t���͕s��)                          *
*      �͂��߂ɂ��̃T�u���[�`���֓n�����X�ł�����                  *
* X1...�e�����ւ̑|�����X(�ԍ��t���͎������ɈقȂ�)                *
* XO...�e�����ւ̑|���O��X(�ԍ��t���͕s��)                          *
*********************************************************************
      SUBROUTINE LBLBND (AT1,AT2,AT3,AT4,AT5,BX,XN,NX,NY,NE,
     $                   X1,XO,A,B,C,D,U,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE)
      DIMENSION XN(NE),X1(NE),BX(NE),XO(NE)
      DIMENSION A(NE),B(NE),C(NE),D(NE),U(NE),Y(NE)
* �ő�J��Ԃ���(NITR)�Ǝ�������l(EITR)�̐ݒ�
      NITR = 200
      EITR = 1.0D-9
* �ɘa�W���̐ݒ� ( OMG=1.0D0�Ȃ�Gauss-Seidel�@ )
      WRITE (6,*) ' Input OMG ( if Gauss-Seidel, OMG=1.0D0 ) ---> '
      READ (5,*) OMG
* B��2��m�����̌v�Z
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + BX(I)**2
   10 CONTINUE
*ITR : �J��Ԃ��񐔂̃J�E���^
      DO 20 ITR=1,NITR
*       x �������ւ̑|�� : �g�[�}�X�@�ɂ��
        INX = 1
        DO 100 IY = 1,NY
          DO 110 IX = 1,NX
            INY = IY + (IX-1)*NY
*           �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X
            A(INX) = AT1(INY)
            B(INX) = AT3(INY)
            C(INX) = AT5(INY)
            D(INX) = BX(INY)
            IF (INY-1.GE.1) THEN
              D(INX)=D(INX)-AT2(INY)*XN(INY-1)
            END IF
            IF (INY+1.LE.NE) THEN
              D(INX)=D(INX)-AT4(INY)*XN(INY+1)
            END IF
*           �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ�
            XO(INY) = XN(INY)
            INX = INX + 1
  110     CONTINUE
  100   CONTINUE
*       Ly=b ������
        U(1) = C(1) / B(1)
        DO 120 J = 2,NE-1
        U(J) = C(J) / ( B(J)-A(J)*U(J-1) )
  120   CONTINUE
        Y(1) = D(1) / B(1)
        DO 130 J = 2,NE
          Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  130   CONTINUE
*       Ux=y ������
        X1(NE) = Y(NE)
        DO 140 J = NE-1,1,-1
          X1(J) = Y(J) - U(J)*X1(J+1)
  140   CONTINUE
        INX = 1
        DO 150 IY = 1,NY
          DO 160 IX = 1,NX
            INY = IY + (IX-1)*NY
*           ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa
            XN(INY)=(1.0D0-OMG)*XO(INY)+OMG*X1(INX)
            INX = INX + 1
  160     CONTINUE
  150   CONTINUE
*       y �������ւ̑|�� : �g�[�}�X�@�ɂ��
        INY = 1
        DO 200 IX = 1,NX
          DO 210 IY = 1,NY
            INX = IX + (IY-1)*NX
*           �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X
            A(INY) = AT2(INY)
            B(INY) = AT3(INY)
            C(INY) = AT4(INY)
            D(INY) = BX(INY)
            IF (INY-NY.GE.1) THEN
              D(INY)=D(INY)-AT1(INY)*XN(INY-NY)
            END IF
            IF (INY+NY.LE.NE) THEN
              D(INY)=D(INY)-AT5(INY)*XN(INY+NY)
            END IF
*           �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ�
            XO(INY) = XN(INY)
            INY = INY + 1
  210     CONTINUE
  200   CONTINUE
*       Ly=b ������
        U(1) = C(1) / B(1)
        DO 220 J = 2,NE-1
          U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  220   CONTINUE
        Y(1) = D(1) / B(1)
        DO 230 J = 2,NE
          Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  230   CONTINUE
*       Ux=y ������
        X1(NE) = Y(NE)
        DO 240 J = NE-1,1,-1
          X1(J) = Y(J) - U(J)*X1(J+1)
  240   CONTINUE
        INY = 1
        DO 250 IX = 1,NX
          DO 260 IY = 1,NY
*           ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa
            XN(INY)=(1.0D0-OMG)*XO(INY)+OMG*X1(INY)
            INY = INY + 1
  260     CONTINUE
  250   CONTINUE
*       RNORM : �c���x�N�g��(R=AX-B)��2��m���� : �T�u���[�`��SORBND�̏ꍇ�������Q��
        RNORM= 0.0D0
        DO 310 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=AT3(I)*XN(I)+AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
          SUM=AT2(I)*XN(I-1)+AT3(I)*XN(I)+AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NE-NY) THEN
            SUM=AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)+AT3(I)*XN(I)
     $         +AT4(I)*XN(I+1)+AT5(I)*XN(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
          SUM=AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)+AT3(I)*XN(I)+AT4(I)*XN(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=AT1(I)*XN(I-NY)+AT2(I)*XN(I-1)+AT3(I)*XN(I)
          END IF
          RNORM= RNORM + (BX(I)-SUM)**2
  310   CONTINUE
*       �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ��������
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) 'Converged : Total ITR = ',ITR
          GO TO 900
        END IF
   20 CONTINUE
*     NITR�܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not converged ! '
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  900 CONTINUE
* ���ʏo��
      DO 320 I = 1,NE
        WRITE(6,2000) I,XN(I)
 2000   FORMAT(' X_9_( ',I3,' ) = ',1PD13.5)
  320 CONTINUE
      RETURN
      END
*********************************************************************
*              �N�����t������Ԗ@ : ��@ 10 - 13
*********************************************************************
*    �����c��(Conjugate Residual)�@�ɂ���Ώ̍s�� A ���܂�        *
*  ���`�V�X�e����@�T�u���[�`�� AX=B                                *
*�m�ϐ��̐����n                                                     *
*    NX:x�����i�q������; NY:y�����i�q������; NE:���i�q�_��=NX*NY    *
*    NITR : �ő唽����; EPSP : ��������l                         *
*�m�z��̐����n                                                     *
*    A(NE,NE) : AX=B �̃t���}�g���b�N�X A                           *
*    B(NE) : ���`�V�X�e�� AX=B �̃x�N�g�� B                         *
*    X(NE) : ���`�V�X�e�� AX=B �̃x�N�g�� X ---> ��������߂�       *
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + ��_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P , AR(NE) : A * R                                *
*********************************************************************
      SUBROUTINE KCR (A,B,X,NE,
     $                R,P,AP,AR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* ��Ɨp�z��
      DIMENSION R(NE),P(NE),AP(NE),AR(NE)
* �ő�J��Ԃ��񐔂̐ݒ�
      NITR = 200
* ��������l(�ϐ� ZANSA �����̒l�ȉ��ɂȂ�� NITR �ȉ��ł������Ɣ��肷��)
      EITR = 1.0D-9
* �z��̃[���N���A
      DO 10 J=1, NE
        R(J) = 0.0D0
        P(J) = 0.0D0
        AP(J) = 0.0D0
        AR(J) = 0.0D0
   10 CONTINUE
* R �� AX ����
      CALL PROFMV (A,X,R,NE)
* B�̃m�����̌v�Z
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* r_{0}��p_{0}(�����l)�̐ݒ�
      DO 30 I = 1, NE
*       r_{0} = B - AX �̌v�Z
        R(I) = B(I) - R(I)
*       p_{0} = r_{0}
        P(I) = R(I)
   30 CONTINUE
* AP�� A p_{0} ����(�ȍ~��AP�� ��(A) �ŋ��߂�)
      CALL PROFMV (A,P,AP,NE)
* �J��Ԃ��v�Z
      DO 40 K = 1,NITR
*       ( r_{k}, A p_{k} )�̌v�Z => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )�̌v�Z => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       ��_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
        ALP = RAP / APAP
        DO 50 I = 1,NE
*         x_{k+1}=x_{k}+��_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-��_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
   50   CONTINUE
*       RNORM : �c���x�N�g��(R=AX-B)��2��m����
        RNORM = 0.0D0
        DO 60 I = 1,NE
          SUM = 0.0D0
          DO 70 J = 1,NE
            SUM = SUM + A(I,J)*X(J)
   70     CONTINUE
          RNORM = RNORM + (B(I)-SUM)**2
   60   CONTINUE
*       �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA �� EITR �ȉ��Ȃ�����Ƃ݂Ȃ��Č��ʂ�\������ 900 ��
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',K
          GO TO 900
        ELSE
*         A r_{k+1} �̌v�Z => AR(NE)
          CALL PROFMV (A,R,AR,NE)
*         ( A r_{k+1}, A p_{k} )�̌v�Z => ARAP
          CALL PROVV( AR,AP,ARAP,NE)
*         ��_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
          BETA = - ARAP / APAP
          DO 80 I = 1, NE
*           p_{k+1} = r_{k+1} + ��_{k} p_{k}
            P(I) = R(I) + BETA*P(I)
*           A p_{k+1} = A r_{k+1} + ��_{k}A p_{k}<-------��(A)
            AP(I) = AR(I) + BETA*AP(I)
   80     CONTINUE
        END IF
   40 CONTINUE
*     NITR�܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not Converged ! '
*     �����Ɣ��肳�ꂽ�Ƃ��̕���_
  900 CONTINUE
*���ʏo��
      DO 90 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_10_( ',I3,' ) = ',1PD13.5)
   90 CONTINUE
      RETURN
      END
*********************************************************************
*    �t���}�g���b�N�X A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`�� AB=C    *
*    NE : ���i�q�_��                                                *
*********************************************************************
      SUBROUTINE PROFMV (A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),C(NE)
      DO 10 I = 1, NE
        S = 0.0D0
        DO 20 J = 1, NE
          S = S + A(I,J)*B(J)
   20   CONTINUE
        C(I) = S
   10 CONTINUE
      RETURN
      END
*********************************************************************
*    �x�N�g�� A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��                 *
*                        AB=C                                       *
*�m�ϐ��̐����n                                                     *
*    NE : ���i�q�_��(�x�N�g�� A,B �̃T�C�Y)                         *
*    C  : A �� B �̐�(�X�J���[)                                     *
*********************************************************************
      SUBROUTINE PROVV(A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(NE),B(NE)
      C = 0.0D0
      DO 10 I=1,NE
        C = C + A(I)*B(I)
   10 CONTINUE
      RETURN
      END
*********************************************************************
*    �����c��(Conjugate Residual)�@�ɂ���Ώ̍s�� A ���܂�        *
*  ���`�V�X�e����@�T�u���[�`��   (�o���h�}�g���b�N�X�p)            *
*                        AX=B                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*  B : ���m�x�N�g��                                                 *
*  X : ���m�x�N�g�� ---> ��������߂�                               *
*�m�ϐ��̐����n                                                     *
*    NX:x�����i�q������; NY:y�����i�q������; NE:���i�q�_��=NX*NY    *
*    NITR : �ő唽����; EPSP : ��������l                         *
*�m�z��̐����n                                                     *
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + ��_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * R                                                 *
*********************************************************************
      SUBROUTINE KCRBND (A1,A2,A3,A4,A5,B,X,NY,NE,
     $                   R,P,AP,AR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),B(NE),X(NE)
* ��Ɨp�z��
      DIMENSION R(NE),P(NE),AP(NE),AR(NE)
* �ő�J��Ԃ��񐔂̐ݒ�
      NITR = 200
* �����������(�ϐ� ZANSA �����̒l�ȉ��ɂȂ�� NITR �ȉ��ł������Ɣ��肷��)
      EITR = 1.0D-9
* �z��̃[���N���A
      DO 10 J=1, NE
        R(J) = 0.0D0
        P(J) = 0.0D0
        AP(J) = 0.0D0
        AR(J) = 0.0D0
   10 CONTINUE
* R �� AX ����
      CALL PROBMV (A1,A2,A3,A4,A5,X,R,NY,NE)
* B�̃m�����̌v�Z
      BNORM = 0.0D0
      DO 20 I = 1,NE
        BNORM = BNORM + B(I)**2
   20 CONTINUE
* r_{0}��p_{0}(�����l)�̐ݒ�
      DO 30 I = 1, NE
*       r_{0} = B - AX �̌v�Z
        R(I) = B(I) - R(I)
*       p_{0} = r_{0}
        P(I) = R(I)
   30 CONTINUE
* AP�� A p_{0} ����(�ȍ~��AP�� ��(A) �ŋ��߂�)
      CALL PROBMV (A1,A2,A3,A4,A5,P,AP,NY,NE)
* �J��Ԃ��v�Z
      DO 40 K = 1,NITR
*       ( r_{k}, A p_{k} )�̌v�Z => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )�̌v�Z => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       ��_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
        ALP = RAP / APAP
        DO 50 I = 1,NE
*         x_{k+1}=x_{k}+��_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-��_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
   50   CONTINUE
*       RNORM : �c���x�N�g��(R=AX-B)��2��m���� : �T�u���[�`��SORBND�̏ꍇ�������Q��
        RNORM= 0.0D0
        DO 60 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NE-NY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   60   CONTINUE
*       �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA �� EITR �ȉ��Ȃ�����Ƃ݂Ȃ��� 700 ��
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',K
          GO TO 700
        END IF
*       ���������̏ꍇ
*       A r_{k+1} �̌v�Z => AR(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,R,AR,NY,NE)
*       ( A r_{k+1}, A p_{k} )�̌v�Z => ARAP
        CALL PROVV( AR,AP,ARAP,NE)
*       ��_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
        BETA = - ARAP / APAP
        DO 70 I = 1, NE
*         p_{k+1} = r_{k+1} + ��_{k} p_{k}
          P(I) = R(I) + BETA*P(I)
*         A p_{k+1} = A r_{k+1} + ��_{k}A p_{k}<-------��(A)
          AP(I) = AR(I) + BETA*AP(I)
   70   CONTINUE
   40 CONTINUE
* NITR �܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not Converged ! '
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  700 CONTINUE
* ���ʏo��
      DO 80 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_11_( ',I3,' ) = ',1PD13.5)
   80 CONTINUE
      RETURN
      END
*********************************************************************
*    �o���h�}�g���b�N�X A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��       *
*                        AB=C                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*�m�ϐ��̐����n                                                     *
*    NY : y�����i�q������                                           *
*                 ---> �o���h�}�g���b�N�X�ƃx�N�g���̐ς̌v�Z�Ŏg�p *
*    NE : ���i�q�_��                                                *
*    C  : A �� B �̐�(�x�N�g��)                                     *
*********************************************************************
      SUBROUTINE PROBMV (A1,A2,A3,A4,A5,B,C,NY,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),B(NE),C(NE)
* �T�u���[�`��SORBND�̏ꍇ�������Q��
      I=1
        C(I) = A3(I)*B(I)+A4(I)*B(I+1)+A5(I)*B(I+NY)
      DO 10 I=2,NY
        C(I) = A2(I)*B(I-1)+A3(I)*B(I)+A4(I)*B(I+1)+A5(I)*B(I+NY)
   10 CONTINUE
      DO 20 I=NY+1,NE-NY
        C(I) = A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)
     $      +A4(I)*B(I+1)+A5(I)*B(I+NY)
   20 CONTINUE
      DO 30 I=NE-NY+1,NE-1
        C(I) = A1(I)*B(I-NY)+A2(I)*B(I-1)+A3(I)*B(I)+A4(I)*B(I+1)
   30 CONTINUE
      I=NE
        C(I) = A1(I)*B(I-NY)+A2(I)*B(I-1 )+A3(I)*B(I)
      RETURN
      END
*********************************************************************
*  Bi-CGSTAB �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*                        AX=B                                       *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : �ő唽����                                            *
*    EPSP : �����������                                            *
*�m�z��̐����n                                                     *
*    T(NE) : t_{k} = r_{k} - ��_{k} A p_{k}                         *
*    X(NE) : x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{K} t_{k}          *
*    R(NE) : r_{k+1} = t_{k} - ��_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P(NE) : p_{k+1} = r_{k+1} + ��_{k} ( p_{k}-��_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * T                                                 *
*********************************************************************
      SUBROUTINE KBICG (A,B,X,NE,
     $                  R,AP,AT,P,S,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NE,NE),B(NE),X(NE)
* ��Ɨp�z��
      DIMENSION R(NE),AP(NE),AT(NE),P(NE),S(NE),T(NE)
* �ő�J��Ԃ��񐔂̐ݒ�
      NITR = 200
* �����������(�ϐ� ZANSA �����̒l�ȉ��ɂȂ�� NITR �ȉ��ł������Ɣ��肷��)
      EITR = 1.0D-9
* B�̃m�����̌v�Z
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* �z��̃[���N���A
      DO 20 J=1, NE
        R(J) = 0.0D0
        AP(J) = 0.0D0
        AT(J) = 0.0D0
        P(J) = 0.0D0
        S(J) = 0.0D0
        T(J) = 0.0D0
   20 CONTINUE
* R �� AX ����
      CALL PROFMV (A,X,R,NE)
*r_{0}��p_{0}(�����l)�C������ s=r_{0} �̐ݒ�
      DO 30 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   30 CONTINUE
* �J��Ԃ��v�Z
      DO 40 J =1,NITR
*       ( s, r_{k} ) �̌v�Z => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} �̌v�Z => AP(NE)
        CALL PROFMV (A,P,AP,NE)
*       ( s, A p_{k} ) �̌v�Z => SAP
        CALL PROVV (S,AP,SAP,NE)
*       ��_{k} = ( s, r_{k} ) / ( s, A p_{k} )
        ALPHA = SR1/SAP
          DO 50 I=1,NE
*           t_{k} = r_{k} - ��_{k} A p_{k}
            T(I) = R(I) - ALPHA*AP(I)
   50     CONTINUE
*       A t_{k} �̌v�Z => AT(NE)
        CALL PROFMV (A,T,AT,NE)
*       ( A t_{k}, t_{k} ) �̌v�Z => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) �̌v�Z => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ��_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
        XI = ATT/ATAT
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ��_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
   60   CONTINUE
* RNORM : �c���x�N�g��(R=AX-B)��2��m����
      RNORM = 0.0D0
      DO 70 I = 1,NE
        SUM = 0.0D0
        DO 80 M = 1,NE
          SUM = SUM + A(I,M)*X(M)
   80   CONTINUE
        RNORM = RNORM + (B(I)-SUM)**2
   70 CONTINUE
*     �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA �� EITR �ȉ��Ȃ�����Ƃ݂Ȃ��� 900 ��
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',J
          GO TO 900
        END IF
*       ���������̏ꍇ : ��_{k}�� p_{k+1} �����߂ČJ��Ԃ��v�Z
*       ( s, r_{k+1} ) �̌v�Z => SR2
        CALL PROVV (S,R,SR2,NE)
*       ��_{k} = ( ��_{k} / ��_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 90 I=1,NE
*         p_{k+1} = r_{k+1} + ��_{k}( p_{k}-��_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   90   CONTINUE
   40 CONTINUE
*     NITR �܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not Converged ! '
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  900 CONTINUE
* ���ʏo��
      DO 100 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_12_( ',I3,' ) = ',1PD13.5)
  100 CONTINUE
      RETURN
      END
*********************************************************************
*  Bi-CGSTAB �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*  (�o���h�}�g���b�N�X�p)                                           *
*                        AX=B                                       *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE)                *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in2d.mac)�ɂĐݒ�                *
*�m�z��̐����n                                                     *
*    T(NE) : t_{k} = r_{k} - ��_{k} A p_{k}                         *
*    X(NE) : x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{K} t_{k}          *
*    R(NE) : r_{k+1} = t_{k} - ��_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P(NE) : p_{k+1} = r_{k+1} + ��_{k} ( p_{k}-��_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP(NE) : A * P                                                 *
*    AT(NE) : A * T                                                 *
*********************************************************************
      SUBROUTINE KBIBND (A1,A2,A3,A4,A5,B,X,NY,NE,
     $                   R,AP,AT,P,S,T)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),B(NE),X(NE)
* ��Ɨp�z��
      DIMENSION R(NE),AP(NE),AT(NE),P(NE),S(NE),T(NE)
* �ő�J��Ԃ��񐔂̐ݒ�
      NITR = 200
* �����������(�ϐ� ZANSA �����̒l�ȉ��ɂȂ�� NITR �ȉ��ł������Ɣ��肷��)
      EITR = 1.0D-9
* B�̃m�����̌v�Z
      BNORM = 0.0D0
      DO 10 I = 1,NE
        BNORM = BNORM + B(I)**2
   10 CONTINUE
* �z��̃[���N���A
      DO 20 J=1, NE
        R(J) = 0.0D0
        AP(J) = 0.0D0
        AT(J) = 0.0D0
        P(J) = 0.0D0
        S(J) = 0.0D0
        T(J) = 0.0D0
   20 CONTINUE
* R �� AX ����
      CALL PROBMV (A1,A2,A3,A4,A5,X,R,NY,NE)
* r_{0}��p_{0}(�����l)�C������ s=r_{0} �̐ݒ�
      DO 30 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   30 CONTINUE
* �J��Ԃ��v�Z
      DO 40 J =1,NITR
*       ( s, r_{k} ) �̌v�Z => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} �̌v�Z => AP(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,P,AP,NY,NE)
*       ( s, A p_{k} ) �̌v�Z => SAP
        CALL PROVV (S,AP,SAP,NE)
*       ��_{k} = ( s, r_{k} ) / ( s, A p_{k} )
        ALPHA = SR1/SAP
        DO 50 I=1,NE
*         t_{k} = r_{k} - ��_{k} A p_{k}
          T(I) = R(I) - ALPHA*AP(I)
   50   CONTINUE
*       A t_{k} �̌v�Z => AT(NE)
        CALL PROBMV (A1,A2,A3,A4,A5,T,AT,NY,NE)
*       ( A t_{k}, t_{k} ) �̌v�Z => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) �̌v�Z => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ��_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
        XI = ATT/ATAT
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ��_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
   60   CONTINUE
*       RNORM : �c���x�N�g��(R=AX-B)��2��m���� : �T�u���[�`��SORBND�̏ꍇ�������Q��
        RNORM= 0.0D0
        DO 70 I = 1,NE
          SUM = 0.0D0
          IF (I.EQ.1) THEN
            SUM=A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.2. AND. I.LE.NY) THEN
            SUM=A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NY+1. AND. I.LE.NE-NY) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
     $         +A4(I)*X(I+1)+A5(I)*X(I+NY)
          ELSE IF (I.GE.NE-NY+1. AND. I.LE.NE-1) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)+A4(I)*X(I+1)
          ELSE IF (I.EQ.NE) THEN
            SUM=A1(I)*X(I-NY)+A2(I)*X(I-1)+A3(I)*X(I)
          END IF
          RNORM= RNORM + (B(I)-SUM)**2
   70   CONTINUE
*       �c���̌v�Z ZANSA = || R || / || B ||
        ZANSA = DSQRT(RNORM/BNORM)
*       ZANSA �� EITR �ȉ��Ȃ�����Ƃ݂Ȃ��� 900 ��
        IF (ZANSA.LE.EITR) THEN
          WRITE (6,*) ' Converged : Total ITR = ',J
          GO TO 900
        END IF
*       ���������̏ꍇ
*       ( s, r_{k+1} ) �̌v�Z => SR2
        CALL PROVV (S,R,SR2,NE)
*       ��_{k} = ( ��_{k} / ��_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 90 I=1,NE
*         p_{k+1} = r_{k+1} + ��_{k}( p_{k}-��_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   90   CONTINUE
   40 CONTINUE
* NITR �܂Ōv�Z���Ă���������
      WRITE (6,*) ' Not Converged ! '
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  900 CONTINUE
* ���ʏo��
      DO 100 I = 1,NE
        WRITE(6,2000) I,X(I)
 2000   FORMAT(' X_13_( ',I3,' ) = ',1PD13.5)
  100 CONTINUE
      RETURN
      END
