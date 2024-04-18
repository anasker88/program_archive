*********************************************************************
* �t�@�C����  �FHSMAC2D.FOR                                         *
* �^�C�g��    �FHSMAC�@�ɂ��2�����M������̓v���O����              *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� ���p���w��                      *
* �����      �F2003.12.25                                          *
* ����        �FFORTRAN (FORTRAN77�ł����s�\)                     *
*********************************************************************
*                                                                   *
*    �{�v���O�����ł́C�Η����͔�ۑ��`��1�����x���㍷�����C�g�U��  *
*  ��2�����x���S������p���ė��U�����Ă���D����ɍ����x�̋ߎ����s��*
*  �ꍇ�́C�K�X�ύX�̂��ƁD                                         *
*  �i�q��������ύX����Ƃ��́CPARAMETER������NX0,NY0�����ׂĕύX�D *
*                                                                   *
*  ���ϐ�                                                           *
*         ��                                                        *
*    ���x �u=(U,V,W),   ���� P,   ���x T                            *
*                                                                   *
*  ����b�������ɂ���                                             *
*        ��                                                         *
*    ���E�u ���O                                                    *
*      ��                                                           *
*    �c�u                    �Q��                                   *
*    �|�|���|���o�{�u�h�r�i���@�u�j�{�a�t�n�~�s(0,1)                *
*    �c��          ~~~~~~            ~~~~~~                         *
*    �c�s            �Q                                             *
*    �|�|���`�k�o�i��  �s�j                                         *
*    �c��  ~~~~~~                                                   *
*    �������ɉ����āC"BIS","BUO","ALP"���`���ė^����D            *
*                                                                   *
*  ���i�q�����ɂ��� (NX=3, NY=3�̗�)                              *
*       ���z�Z�� �������E                  �E�����E ���z�Z��        *
*             ��     ��                         �� ��               *
*           +--------�{-------+--------+--------�{-------+          *
*           |P(0,NY+1)        |        |        �aP(NX+1,NY+1)      *
* ���z�Z����|   �E   ��  �E   |   �E   |   �E   ��  �E   |�����z�Z��*
*           |      U(0,NY+1)  |        |      U(NX,NY+1) |          *
* ������E��+===��===�{=======+========+========�{==��===+��������E*
* (IY=NY)   |V(0,NY) �a       |        |        �aV(NX+1,NY)        *
*           |   �E   �a  �E   |   �E   |   �E   �a  �E   |          *
*           |        �a       |        |        �a       |          *
*           +--------�{-------+--------+--------�{-------+          *
*           |        �a       |        |        �a       |          *
*           |   �E   �a  �E   |   �E   |   �E   �a  �E   |          *
*           |        �a       |        |        �a       |          *
*           +--------�{-------+--------+--------�{-------+          *
*           |        �a       |        |        �a       |          *
*           |   �E   �a  �E   |   �E   |   �E   �a  �E   |          *
*           |  V(0,0)�a       |        |        �aV(NX+1,0)         *
* �������E��+===��===�{=======+========+========�{==��===+���������E*
* (IY=0)    | P(0,0) �a       |        |        �aP(NX+1,0)         *
* ���z�Z����|   �E   ��  �E   |   �E   |   �E   ��  �E   |�����z�Z��*
*           |      U(0,0)     |        |      U(NX,0)    |          *
*     ��    +--------�{-------+--------+--------�{-------+          *
*     ��       ��    ��                         ��    ��            *
*  �� �b  ���z�Z��  �������E                  �E�����E ���z�Z��     *
*  �� �b            (IX=0)                    (IX=NX)               *
*     �{�|����                                                      *
*                                                                   *
*  ���X�^�b�K�[�h���b�V���ɂ���                                   *
*                               ��    DX     ��                     *
*    +------------+-------------+-------------+                     *
*    |            |             |             | ��                  *
*    |            |    P(i,j+1) |             |                     *
*    |            |      �E     |             | DY                  *
*    |            |             |             |                     *
*    |            |    V(i,j)   |             | ��                  *
*    +------------+------��-----+-------------+                     *
*    |            |             |             |                     *
*    |   P(i-1,j) |    P(i,j)   |   P(i+1,j)  | �� U ��`�_         *
*    |     �E     ��     �E     ��     �E     | �� V ��`�_         *
*    |         U(i-1,j)       U(i,j)          | �E P,T ��`�_       *
*    |            |    V(i,j-1) |             | (���j               *
*    +------------+------��-----+-------------+ �v���O��������T�́C *
*    |            |             |             | �{�����ł̓��ƂȂ���*
*    |            |    P(i,j-1) |             | ����D              *
*    |            |      �E     |             |                     *
*    |            |             |             |                     *
*    |            |             |             |                     *
*    +------------+-------------+-------------+                     *
*                                                                   *
*  ���p�����[�^�[�t�@�C���i�����Q���D�������j�ɂ���               *
*    �v���O���������s����ƁC�h�����Q���D�������h�Ƃ����p�����[�^   *
*    �t�@�C����ǂ݂ɍs���̂ŁC���炩���ߍ쐬���Ă����D             *
*                                                                   *
*  "in2d.mac"�̃��X�g                                               *
* 1: U.NEW ....�t�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 2: V.NEW ....�u�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 3: P.NEW ....�o�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 4: T.NEW ....�s�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 5: U.OLD ....�p���v�Z�̂t�̓��̓f�[�^                             *
* 6: V.OLD ....�p���v�Z�̂u�̓��̓f�[�^                             *
* 7: P.OLD ....�p���v�Z�̂o�̓��̓f�[�^                             *
* 8: T.OLD ....�p���v�Z�̂s�̓��̓f�[�^                             *
* 9: UVT.NEW ....Tecplot�p�f�[�^                                    *
*10: +============================+                                 *
*11: | ITYPE   1===>   isothermal |                                 *
*12: |         2===>nonisothermal |                                 *
*13: +============================+                                 *
*14: -----------------------------------------------------------    *
*15: ITYPE     ICYCLE   NITR    NCYCLE                              *
*16: 2         0        10000   1000                                *
*17: -----------------------------------------------------------    *
*18: EPSP      OMG                                                  *
*19: 1.0e-3    1.7e+0                                               *
*20: -----------------------------------------------------------    *
*21: DT        RE       PR       GR                                 *
*22: 1.0e-4    0.0e+0   7.1e-1   1.0e+5                             *
*23: -----------------------------------------------------------    *
*24: DLX       DLY      IRELP   METHOD                              *
*25: 1.0e+0    1.0e+0   0       5                                   *
*                                                                   *
*    ITYPE.....�P�P�C�P�Q�s���Q��                                   *
*    ICYCLE.....�v�Z�J�n�̃T�C�N�����i����t=ICYCLE*DT�j             *
*    NITR.....���͌v�Z�̂��߂́C�P�T�C�N��������̍ő唽����      *
*    NCYCLE.....�v�Z�I���T�C�N����                                  *
*                             ��                                    *
*    EPSP.....��������l�i���E�u ��EPSP�𖞑�����܂ŁC�����v�Z�ɂ� *
*             ���āC���͏���v�Z����D�j                            *
*    OMG.....���͌v�Z�̂��߂̊ɘa�W�� , DT.....���ԍ���             *
*    RE.....���C�m���Y�� , PR.....�v�����g���� , GR.....�O���X�z�t��*
*    DLX.....��͗̈�̉���(DLX=NX*DX)                              *
*    DLY.....��͗̈�̍���(DLY=NY*DY)                              *
*            �i�q���i���Ԋu�jDX,DY�̓v���O�����̒��ŋ��߂�D        *
*    IRELP...���͂̊�l�̐ݒ�(0:�s��Ȃ�; 1:�s��)                 *
*    METHOD..SMAC�@�ɂ����Ă̂ݗL��                                 *
*    �p�����[�^�t�@�C���̐��l�ɂ��āC                             *
*    FORTRAN �v���O����.....�ł���Δ{���x�����ŗ^����"1.0D0,1.0d0" *
*            C�v���O�����Ƌ��p������"1.0E0,1.0e0"�Ƃ��Ă����͂Ȃ� *
*    C �v���O����....."1.0E0,1.0e0"�Ƃ��ė^����                     *
*    TECPLOT�p�f�[�^�����������o�̓t�@�C���͏����Ȃ��`���ŁC�g�p    *
*    ����R���p�C���[�Ɉˑ�����D�R���p�C���[�Ɉˑ����Ȃ��`���ɂ��� *
*    �ɂ́C�e�ʂ͑����邪�����t���`���ɕύX����΂悢�D             *
*                                                                   *
*  �����͂̑��ΐ��ɂ���                                           *
*    ���͂̑��ΐ��ݒ�̕ϐ���:IRELP                                 *
*    0:���͂̊��݂��Ȃ��D                                       *
*    1:���͂̊��݂���D                                         *
*      �{�v���O�����ł�PO(1,1)=0�ƂȂ�悤�ɐݒ肵�Ă���D          *
*      �T�u���[�`�� PRESS ���Q�ƁD                                  *
*                                                                   *
*  ���ϐ��E�z��̐���                                               *
*    ICYCLE -----> ���Ԑi�s�̂��߂̃J�E���^                         *
*    ITR    -----> ���͌v�Z�̂��߂̔����񐔂̃J�E���^               *
*    IX,IY  -----> ��̐}���Q��                                     *
*    UO,VO,TO----->���͂̔����v�Z���s���O�̒l                       *
*    UN,VN,TN----->���������V���Ȉ��͂�p���Čv�Z���ꂽ�l           *
*                                                                   *
*********************************************************************
      PROGRAM HSMAC2D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      CHARACTER FNAME(10)*20
*x�����̊i�q������
      NX  = NX0
*y�����̊i�q������
      NY  = NY0
*�p�����[�^�t�@�C���̃I�[�v��
      OPEN (10,FILE='IN2D.MAC',STATUS='OLD')
*�o�̓t�@�C�����̓ǂݍ���
      DO 10 I = 1,9
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*U�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
*�����Ȃ��`���̓R���p�C���[�Ɉˑ�����̂Œ���
      OPEN (11, FILE=FNAME(1), STATUS='NEW', FORM='UNFORMATTED')
*V�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (12, FILE=FNAME(2), STATUS='NEW', FORM='UNFORMATTED')
*P�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (13, FILE=FNAME(3), STATUS='NEW', FORM='UNFORMATTED')
*T�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (14, FILE=FNAME(4), STATUS='NEW', FORM='UNFORMATTED')
*in2d.mac���̃R�����g�s(10-15�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) ITYPE, ICYCLE, NITR, NCYCLE
*in2d.mac���̃R�����g�s(17-18�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) EPSP, OMG
*in2d.mac���̃R�����g�s(20-21�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DT,RE, PR, GR
*in2d.mac���̃R�����g�s(23-24�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DLX, DLY, IRELP, METHOD
*�p���̌v�Z�̏ꍇ
      IF (ICYCLE.NE.0) THEN
*U�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (15, FILE=FNAME(5), STATUS='OLD', FORM='UNFORMATTED')
*V�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (16, FILE=FNAME(6), STATUS='OLD', FORM='UNFORMATTED')
*P�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (17, FILE=FNAME(7), STATUS='OLD', FORM='UNFORMATTED')
*T�f�[�^�t�@�C���̃I�[�v��(������ł�T=0.0�̃f�[�^��ǂݍ���)(�����Ȃ��`��)
        OPEN (18, FILE=FNAME(8), STATUS='OLD', FORM='UNFORMATTED')
      END IF
*x�����̊i�q��
      DX  = DLX / FLOAT(NX)
*y�����̊i�q��
      DY  = DLY / FLOAT(NY)
*�^�����������̊g�U���̌W��(�����ł�Pr)
      VIS = PR
*�G�l���M�[���������̊g�U���̌W��(�����ł�1)
      ALP = 1.0D+0
*���͍��̌W��(�����ł� Gr * Pr**2)
      BUO = GR * PR**2
*������Ȃ畂�͍��̌W���̓[���ɐݒ�
      IF (ITYPE.EQ.1) BUO = 0.0D0
*�����l�̐ݒ�
      CALL CINITI
*���Ԑi�s�̂��߂̖߂�_
  700 CONTINUE
*���Ԑi�s
      CALL ADV
*���x��̌v�Z
      CALL CALVEL
*���͌v�Z�̔����񐔂�1�ɏ�����
      ITR = 1
*���͔����̂��߂̖߂�_
  710 CONTINUE
*���͌v�Z�������������ǂ����̃p�����[�^IFLG��������
*IFLG -> 0:���� 1:���U(�ݒ肳�ꂽ���e��NITR�ȉ��ŉ��������Ȃ�)
      IFLG = 0
*���͏�̌v�Z
      CALL PRESS
*     Newton�@�ɂ�鈳�͏�̌v�Z�����������Ƃ�
      IF ( IFLG. EQ. 0 ) THEN
*       �񓙉���v�Z�̏ꍇ
        IF ( ITYPE.EQ.2 ) THEN
*         ���x����v�Z
          CALL CALTEM
        END IF
*     ���͏�̌v�Z���������Ă��Ȃ��Ƃ�
      ELSE IF ( IFLG. EQ. 1 ) THEN
*       ���͌v�Z�̔����񐔂����炩���ߐݒ肳�ꂽ�ő�lNITR��菬�����Ƃ�
        IF ( ITR. LT. NITR ) THEN
*         ����ɔ������J��Ԃ�
          ITR = ITR + 1
          GO TO 710
*       ���͌v�Z�̔����񐔂�NITR�ƂȂ����甭�U�Ƃ݂Ȃ��Čv�Z�I��
        ELSE
          WRITE (6,*) ' NOT CONVERGE ! '
C         �f�[�^���o�͂��ċ����I��
          CALL PROUT
          GO TO 900
        END IF
      END IF
*     ���Ԑi�s�J�E���^(ICYCLE)��NCYCLE��菬������
      IF ( ICYCLE. LT. NCYCLE ) THEN
        GO TO 700
*     ���Ԑi�s�J�E���^��NCYCLE�ɂȂ�����v�Z�I��
      ELSE
        CALL PROUT
      END IF
  900 CONTINUE
*Tecplot�p�f�[�^�̏o��
      CALL TECPLT(FNAME(9))
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      STOP
      END
*********************************************************************
*                        �����ݒ�
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*�V�K�v�Z�̏ꍇ
      IF ( ICYCLE. EQ. 0 ) THEN
*       U�̏����l�ݒ�
        DO 10 IX = 0,NX
          DO 20 IY = 0,NY+1
            UN(IX,IY) = 0.0D0
   20     CONTINUE
   10   CONTINUE
*       V�̏����l�ݒ�
        DO 30 IX = 0,NX+1
          DO 40 IY = 0,NY
            VN(IX,IY) = 0.0D0
   40     CONTINUE
   30   CONTINUE
*       P�̏����l�ݒ�
        DO 50 IX = 0,NX+1
          DO 60 IY = 0,NY+1
            PO(IX,IY) = 0.0D0
   60     CONTINUE
   50   CONTINUE
*----------------------------------------------------------------------*
*�i���Ӂj���͍��̌v�Z�ŉ��x�̔z����g�p���Ă���̂œ�����ł�T=0�Ƃ��� *
* �������������͐ݒ肷��K�v������D�[���ȊO�̒l������ƕ��͍����v�Z *
* �����\��������̂Œ��ӁD                                         *
*----------------------------------------------------------------------*
*       T�̏����l�ݒ�(�̈���͍���(+0.5)�ƒቷ(-0.5)�̒��ԉ��x)
        DO 61 IX = 0,NX+1
          DO 62 IY = 0,NY+1
            TN(IX,IY) = 0.0D0
   62     CONTINUE
   61   CONTINUE
*       T�̋��E�F�E�ʁi��p�jT=-0.5
        DO 70 IY = 0,NY+1
          TN(NX+1,IY) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY)
   70   CONTINUE
*       T�̋��E�F���ʁi���M�jT=+0.5
        DO 80 IY = 0,NY+1
          TN(0,IY) = 2.0D0 * ( +0.5D0 ) - TN(1,IY)
   80   CONTINUE
*       T�̋��E�F��ʁi�f�M�j
        DO 90 IX = 1,NX
          TN(IX,NY+1) = TN(IX,NY)
   90   CONTINUE
*       T�̋��E�F���ʁi�f�M�j
        DO 95 IX = 1,NX
          TN(IX,0) = TN(IX,1)
   95   CONTINUE
*�p���v�Z�i���łɂ���v�Z���ʂ���X�^�[�g�j�̏ꍇ
      ELSE
*       U�f�[�^�t�@�C������̓ǂݍ���[Unit No.=15](�����Ȃ��`��)
        READ (15) UN
*       V�f�[�^�t�@�C������̓ǂݍ���[Unit No.=16](�����Ȃ��`��)
        READ (16) VN
*       P�f�[�^�t�@�C������̓ǂݍ���[Unit No.=17](�����Ȃ��`��)
        READ (17) PO
*---------------------------------------------------------------------*
*    (����) ������̌v�Z�ł�T(=0)�̃t�@�C����ǂݍ��ޕK�v������       *
*---------------------------------------------------------------------*
*       T�f�[�^�t�@�C������̓ǂݍ���[Unit No.=18](�����Ȃ��`��)
        READ (18) TN
        CLOSE (15)
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
      END IF
      RETURN
      END
*********************************************************************
*                         ���Ԑi�s
*********************************************************************
      SUBROUTINE ADV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      TIME   = DT*FLOAT(ICYCLE)
      ICYCLE = ICYCLE + 1
*���Ԑi�s�J�E���^(ICYCLE)��100�񖈂ɕ\��
      IF (MOD(ICYCLE,100).EQ.0) THEN
        WRITE (6,2000) ICYCLE
 2000   FORMAT ('  CYC = ',I8)
      END IF
*--------------------------------------------------------------------
* UN -> UO : �K�v�Ȃ����ւ���O��UN��UO����ϓ��ʂ����߂�
* UN : �O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����
* UO : �V�������ԃX�e�b�v�ł̏����l�DUN��ۑ��D
*--------------------------------------------------------------------
      DO 70 IX = 0,NX
        DO 80 IY = 0,NY+1
          UO(IX,IY) = UN(IX,IY)
   80   CONTINUE
   70 CONTINUE
*--------------------------------------------------------------------
* VN -> VO : �K�v�Ȃ����ւ���O��VN��VO����ϓ��ʂ����߂�
* VN : �O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����
* VO : �V�������ԃX�e�b�v�ł̏����l�DVN��ۑ��D
*--------------------------------------------------------------------
      DO 90 IX = 0,NX+1
        DO 100 IY = 0,NY
          VO(IX,IY) = VN(IX,IY)
  100   CONTINUE
   90 CONTINUE
*--------------------------------------------------------------------
* TN -> TO : �K�v�Ȃ����ւ���O��TN��TO����ϓ��ʂ����߂�
* TN : �O�̎��ԃX�e�b�v�ł̒l
* TO : �V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D
*--------------------------------------------------------------------
      DO 110 IX = 0,NX+1
        DO 120 IY = 0,NY+1
          TO(IX,IY) = TN(IX,IY)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                    ���x��̌v�Z
*********************************************************************
      SUBROUTINE CALVEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*--------------------------------------------------------------------
*                     U(IX,IY)�̌v�Z
*--------------------------------------------------------------------
      DO 10 IX = 1,NX-1
        DO 20 IY = 1,NY
*       VV��U(IX,IY)�ɂ�����V�̕�Ԓl
        VV = (  VO(IX,IY  )+VO(IX+1,IY  )
     $         +VO(IX,IY-1)+VO(IX+1,IY-1) )/4.0D0
*       �Η���(CNVUX,CNVUY)���P�����x���㍷���ɂČv�Z
        IF ( UO(IX,IY). GE. 0.0D0 ) THEN
          CNVUX = UO(IX,IY)*( UO(IX,IY) - UO(IX-1,IY) ) / DX
        ELSE IF ( UO(IX,IY). LT. 0.0D0 ) THEN
          CNVUX = UO(IX,IY)*( UO(IX+1,IY) - UO(IX,IY) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY) - UO(IX,IY-1) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY+1) - UO(IX,IY) ) / DY
        END IF
*       x�����̕��͍�(BUOU)�̓[��
        TU = 0.0D0
        BUOU = BUO * TU
*       �g�U��(DIFU)�̌v�Z
        DIFU = VIS*(
     $        ( UO(IX-1,IY)-2.0D0*UO(IX,IY)+UO(IX+1,IY) )/DX**2
     $       +( UO(IX,IY-1)-2.0D0*UO(IX,IY)+UO(IX,IY+1) )/DY**2
     $              )
*       ���̑��x(U)�̌v�Z
        UN(IX,IY) = UO(IX,IY)
     $  + DT*( -CNVUX-CNVUY+DIFU+BUOU+( PO(IX,IY)-PO(IX+1,IY) )/DX )
   20   CONTINUE
   10 CONTINUE
*--------------------------------------------------------------------
*                       V(IX,IY)�̌v�Z
*--------------------------------------------------------------------
      DO 30 IX = 1,NX
        DO 40 IY = 1,NY-1
*       UU��V(IX,IY)�ɂ�����U�̕�Ԓl
        UU = (  UO(IX-1,IY  )+UO(IX,IY  )
     $         +UO(IX-1,IY+1)+UO(IX,IY+1) )/4.0D0
*       �Η���(CNVVX,CNVVY)���P�����x���㍷���ɂČv�Z
        IF ( UU. GE. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX,IY) - VO(IX-1,IY) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX+1,IY) - VO(IX,IY) ) / DX
        END IF
        IF ( VO(IX,IY). GE. 0.0D0 ) THEN
          CNVVY = VO(IX,IY)*( VO(IX,IY) - VO(IX,IY-1) ) / DY
        ELSE IF ( VO(IX,IY). LT. 0.0D0 ) THEN
          CNVVY = VO(IX,IY)*( VO(IX,IY+1) - VO(IX,IY) ) / DY
        END IF
*       ���͍�(BUOV)�̌v�Z
        TV = ( TO(IX,IY) + TO(IX,IY+1) )/2.0D0
        BUOV = BUO*TV
*       �g�U��(DIFV)�̌v�Z
        DIFV = VIS*(
     $          ( VO(IX-1,IY)-2.0D0*VO(IX,IY)+VO(IX+1,IY) )/DX**2
     $         +( VO(IX,IY-1)-2.0D0*VO(IX,IY)+VO(IX,IY+1) )/DY**2
     $             )
*       ���̑��x(V)�̌v�Z
        VN(IX,IY) = VO(IX,IY)
     $  + DT*(-CNVVX-CNVVY+DIFV+BUOV+(PO(IX,IY)-PO(IX,IY+1))/DY )
   40   CONTINUE
   30 CONTINUE
*���x�̋��E�����̏���
      CALL VELBND
      RETURN
      END
*********************************************************************
*                        ���͏�̌v�Z
*********************************************************************
      SUBROUTINE PRESS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      IXMAX = 0
      IYMAX = 0
      DMAX = 0.0D0
*P(IX,IY)�̌v�Z
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DEL = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 )
          DIV = ( UN(IX,IY) - UN(IX-1,IY  ) )/DX
     $        + ( VN(IX,IY) - VN(IX  ,IY-1) )/DY
          IF ( DABS(DIV). GE. DABS(DMAX) ) THEN
            IXMAX = IX
            IYMAX = IY
            DMAX = DIV
          END IF
          DELP = - OMG * DIV / DEL
          PO(IX  ,IY  ) = PO(IX  ,IY  ) + DELP
          UN(IX  ,IY  ) = UN(IX  ,IY  ) + DT/DX*DELP
          UN(IX-1,IY  ) = UN(IX-1,IY  ) - DT/DX*DELP
          VN(IX,  IY  ) = VN(IX,  IY  ) + DT/DY*DELP
          VN(IX,  IY-1) = VN(IX,  IY-1) - DT/DY*DELP
   20   CONTINUE
   10 CONTINUE
*
* ���͂̑��ΐ��Ɋւ��鏈��(IRELP=1�Ȃ�ȉ��̏������s��)
      IF (IRELP.EQ.1) THEN
        POSTN = PO(1,1)
        DO 30 IX = 1,NX
          DO 40 IY = 1,NY
            PO(IX,IY) = PO(IX,IY)-POSTN
   40     CONTINUE
   30   CONTINUE
      END IF
*
* IFLG=1�Ȃ�C�A���̎��𖞂����Ă��Ȃ��Ɣ��肵�Ăш��͌v�Z���s���D
      IF ( DABS(DMAX). GE. EPSP ) IFLG = 1
* ���͌v�Z�̉񐔂�100�񂲂Ƃɕ\��
      IF (MOD(ITR,100).EQ.0) THEN
        WRITE (6,2000) ITR, IXMAX, IYMAX, DMAX
 2000   FORMAT (' Iteration=',I8,'   Div(max)(',2I6,')=',1PE13.5)
      END IF
*�V���ɓ���ꂽ���x��p���ċ��E��������������
      CALL VELBND
      RETURN
      END
*********************************************************************
*                      ���x��̌v�Z
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*T(IX,IY)�̌v�Z
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
*         UUT,VVT�͂��ꂼ��T(IX,IY)�ɂ�����U,V�̕�Ԓl
          UUT = ( UO(IX,IY) + UO(IX-1,IY  ) ) / 2.0D0
          VVT = ( VO(IX,IY) + VO(IX  ,IY-1) ) / 2.0D0
*         �Η���(CNVTX,CNVTY)���P�����x���㍷���ɂČv�Z
          IF ( UUT. GE. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX,IY) - TO(IX-1,IY) ) / DX
          ELSE IF ( UUT. LT. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX+1,IY) - TO(IX,IY) ) / DX
          END IF
          IF ( VVT. GE. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY) - TO(IX,IY-1) ) / DY
          ELSE IF ( VVT. LT. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY+1) - TO(IX,IY) ) / DY
          END IF
*         �g�U��(DIFT)�̌v�Z
          DIFT = ALP*(
     $       +( TO(IX-1,IY)-2.0D0*TO(IX,IY)+TO(IX+1,IY) )/DX**2
     $       +( TO(IX,IY-1)-2.0D0*TO(IX,IY)+TO(IX,IY+1) )/DY**2
     $                )
*         ���̎��Ԃ�T�̌v�Z
          TN(IX,IY) = TO(IX,IY) + DT*( -CNVTX-CNVTY+DIFT )
   20   CONTINUE
   10 CONTINUE
*���E�����̏���
      CALL TBND
      RETURN
      END
*********************************************************************
*                    ���x�̋��E�����̏���
*********************************************************************
      SUBROUTINE VELBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*U�i�E�ʁj
      DO 10 IY = 1,NY
        UN(NX,IY) = 0.0D0
   10 CONTINUE
*U�i���ʁj
      DO 20 IY = 1,NY
        UN(0,IY) = 0.0D0
   20 CONTINUE
*U�i��ʁj
      DO 30 IX = 0,NX
        UN(IX,NY+1) = -UN(IX,NY)
   30 CONTINUE
*U�i���ʁj
      DO 40 IX = 0,NX
        UN(IX,0) = -UN(IX,1)
   40 CONTINUE
*V�i�E�ʁj
      DO 50 IY = 1,NY-1
        VN(NX+1,IY) = -VN(NX,IY)
   50 CONTINUE
*V�i���ʁj
      DO 60 IY = 1,NY-1
        VN(0,IY) = -VN(1,IY)
   60 CONTINUE
*V�i��ʁj
      DO 70 IX = 0,NX+1
        VN(IX,NY) = 0.0D0
   70 CONTINUE
*V�i���ʁj
      DO 80 IX = 0,NX+1
        VN(IX,0) = 0.0D0
   80 CONTINUE
      RETURN
      END
*********************************************************************
*                  ���x�̋��E�����̏���
*********************************************************************
      SUBROUTINE TBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
*�E��
      DO 10 IY = 0,NY+1
        TN(NX+1,IY) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY)
   10 CONTINUE
*����
      DO 20 IY = 0,NY+1
        TN(0,IY) = 2.0D0 * ( +0.5D0 ) - TN(1,IY)
   20 CONTINUE
*���
      DO 30 IX = 1,NX
        TN(IX,NY+1) = TN(IX,NY)
   30 CONTINUE
*����
      DO 40 IX = 1,NX
        TN(IX,0) = TN(IX,1)
   40 CONTINUE
      RETURN
      END
*********************************************************************
*                       �f�[�^�o��
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      WRITE (11) UN
      WRITE (12) VN
      WRITE (13) PO
      WRITE (14) TN
      RETURN
      END
*********************************************************************
*                       Tecplot�p�f�[�^�o��
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20 )
      COMMON / D1 / NX,NY
      COMMON / D2 / DX,DY,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 / UO(0:NX0,  0:NY0+1),UN(0:NX0  ,0:NY0+1),
     $                  VO(0:NX0+1,0:NY0  ),VN(0:NX0+1,0:NY0  ),
     $                  PO(0:NX0+1,0:NY0+1),
     $                  TO(0:NX0+1,0:NY0+1),TN(0:NX0+1,0:NY0+1)
      CHARACTER FNAME*20
*
      OPEN (21,FILE=FNAME,STATUS='NEW')
      WRITE (21,*) 'VARIABLES = "X", "Y", "U", "V", "T"'
      NX1 = NX+1
      NY1 = NY+1
      WRITE (21,4000) NX1,NY1
 4000 FORMAT (1H ,'ZONE I=',I3,',J=',I3,', F=POINT')
      DO 10 IY = 0,NY
        DO 20 IX = 0,NX
           X = DX * FLOAT(IX)
           Y = DY * FLOAT(IY)
           U = ( UN(IX  ,IY)+UN(IX  ,IY+1) )/2.0D0
           V = ( VN(IX  ,IY)+VN(IX+1,IY  ) )/2.0D0
           T = ( TN(IX  ,IY)  +TN(IX+1,IY  )
     $          +TN(IX  ,IY+1)+TN(IX+1,IY+1) )/4.0D0
          WRITE (21,4010) X,Y,U,V,T
 4010     FORMAT (1H ,5(1PE11.3))
   20   CONTINUE
   10 CONTINUE
      CLOSE(21)
      RETURN
      END
