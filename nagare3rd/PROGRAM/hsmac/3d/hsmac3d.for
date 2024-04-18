*********************************************************************
* �t�@�C����  �FHSMAC3D.FOR                                         *
* �^�C�g��    �FHSMAC�@�ɂ��3�����M������̓v���O����              *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� ���p���w��                      *
* �����      �F2003.12.25                                          *
* ����        �FFORTRAN (FORTRAN77�ł����s�\)                     *
*********************************************************************
*                                                                   *
*    �{�v���O�����ł́C�Η����͔�ۑ��`��1�����x���㍷�����C�g�U��  *
*  ��2�����x���S������p���ė��U�����Ă���D����ɍ����x�̋ߎ����s��*
*  �ꍇ�́C�K�X�ύX�̂��ƁD                                         *
*  �i�q��������ύX����Ƃ���PARAMETER������NX0,NY0,NZ0�����ׂĕύX.*
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
*    �|�|���|���o�{�u�h�r�i���@�u�j�{�a�t�n�~�s(0,1,0)              *
*    �c��          ~~~~~~            ~~~~~~                         *
*    �c�s            �Q                                             *
*    �|�|���`�k�o�i��  �s�j                                         *
*    �c��  ~~~~~~                                                   *
*    �������ɉ����āC"BIS","BUO","ALP"���`���ė^����D            *
*                                                                   *
*  ���i�q�����ɂ��� (NX=3, NY=3 �̗�)                             *
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
*           +--------�{-------+--------+--------�{-------+          *
*              ��    ��                         ��    ��            *
*         ���z�Z��  �������E                  �E�����E ���z�Z��     *
*                   (IX=0)                    (IX=NX)               *
*                                                                   *
*         Y                                                         *
*         ��                                                        *
*         �b                                                        *
*      �� �b                                                        *
*      �� �b                                                        *
*         �{�|�|��X                                                 *
*        /                                                          *
*       /                                                           *
*      /                                                            *
*     Z                                                             *
*                                                                   *
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
*    |         U(i-1,j)       U(i,j)          | �E P,T ��`         *
*    |            |    V(i,j-1) |             | (���j               *
*    +------------+------��-----+-------------+ �v���O��������T�́C *
*    |            |             |             | �{�����ł̓��ƂȂ���*
*    |            |    P(i,j-1) |             | ����D              *
*    |            |      �E     |             |                     *
*    |            |             |             |                     *
*    |            |             |             |                     *
*    +------------+-------------+-------------+                     *
*                                                                   *
*             Y                                                     *
*             +--------------------------------------- +            *
*            /|         V(IX,IY,IZ)                   /|            *
*           / |             ��                       / |            *
*          /  |             �b                      /  |            *
*         /   |             �b W(IX,IY,IZ-1)       /   |            *
*        /    |             �b /                  /    |            *
*       /     |             �|/                  /     |            *
*      /      |              /                  /      |            *
*     /       |             �|                 /       |            *
*    /        |       P,T(IX,IY,IZ)           /        |            *
*   +----------------------------------------+         |            *
*   |   |---> |            ��                |  |--->U(IX,IY,IZ)    *
*   | U(IX-1,IY,IZ)--------------------------|---------+ X          *
*   |        /0          /                   |        /             *
*   |       /           /   ��               |       /              *
*   |      /           /    �b               |      /               *
*   |     /           �|    �b               |     /                *
*   |    /   W(IX,IY,IZ)    �b               |    /                 *
*   |   /                   �|               |   /                  *
*   |  /                   V(IX,IY-1,IZ)     |  /                   *
*   | /                                      | /                    *
*   |/                                       |/                     *
*   +----------------------------------------+                      *
*  Z                                                                *
*                                                                   *
*                                                                   *
*  ���p�����[�^�[�t�@�C���i�����R���D�������j�ɂ���               *
*    �v���O���������s����ƁC�h�����R���D�������h�Ƃ����p�����[�^   *
*    �t�@�C����ǂ݂ɍs���̂ŁC���炩���ߍ쐬���Ă����D             *
*                                                                   *
*  "in3d.mac"�̃��X�g                                               *
* 1: U.NEW ....�t�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 2: V.NEW ....�u�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 3: W.NEW ....�v�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 4: P.NEW ....�o�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 5: T.NEW ....�s�̌v�Z���ʂ̏o�̓t�@�C����                         *
* 6: U.OLD ....�p���v�Z�̂t�̓��̓f�[�^                             *
* 7: V.OLD ....�p���v�Z�̂u�̓��̓f�[�^                             *
* 8: W.OLD ....�p���v�Z�̂v�̓��̓f�[�^                             *
* 9: P.OLD ....�p���v�Z�̂o�̓��̓f�[�^                             *
*10: T.OLD ....�p���v�Z�̂s�̓��̓f�[�^                             *
*11: UVWT.NEW ....Tecplot�p�f�[�^                                   *
*12: +============================+                                 *
*13: | ITYPE   1===>   isothermal |                                 *
*14: |         2===>nonisothermal |                                 *
*15: +============================+                                 *
*16: -----------------------------------------------------------    *
*17: ITYPE     ICYCLE   NITR    NCYCLE                              *
*18: 2         0        10000   10000                               *
*19: -----------------------------------------------------------    *
*20: EPSP      OMG                                                  *
*21: 1.0e-3    1.7e+0                                               *
*22: -----------------------------------------------------------    *
*23: DT        RE       PR       GR                                 *
*24: 1.0e-4    0.0e+0   7.1e-1   1.0e+5                             *
*25: -----------------------------------------------------------    *
*26: DLX       DLY      DLZ      IRELP   METHOD                     *
*27: 1.0e+0    1.0e+0   1.0e0    0       5                          *
*                                                                   *
*    ITYPE.....�P�R�C�P�S�s���Q��                                   *
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
*    DLZ.....��͗̈�̍���(DLZ=NZ*DZ)                              *
*            �i�q���i���Ԋu�jDX,DY,DZ�̓v���O�����̒��ŋ��߂�D     *
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
*      �{�v���O�����ł�PO(1,1,1)=0�ƂȂ�悤�ɐݒ肵�Ă���D        *
*      �T�u���[�`�� PRESS ���Q�ƁD                                  *
*                                                                   *
*  ���ϐ��E�z��̐���                                               *
*    ICYCLE........�^�C���X�e�b�v�̂��߂̃J�E���^                   *
*    ITR...........���͌v�Z�̂��߂̔����񐔂̃J�E���^               *
*    IX,IY,IZ......�X�J���[�ʂ̒�`�_�̂��C�����W�̓Y��             *
*    UO,VO,W0,TO......���͂̔����v�Z���s���O�̒l                    *
*    UN,VN,WN,TN......���������V���Ȉ��͂�p���Čv�Z���ꂽ���x      *
*                                                                   *
*********************************************************************
      PROGRAM HSMAC3D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      CHARACTER FNAME(11)*20
*x�����̊i�q������
      NX  = NX0
*y�����̊i�q������
      NY  = NY0
*Z�����̊i�q������
      NZ  = NZ0
*�p�����[�^�t�@�C���̃I�[�v��
      OPEN (10,FILE='IN3D.MAC',STATUS='OLD')
*�o�̓t�@�C�����̓ǂݍ���
      DO 10 I = 1,11
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*U�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
*�����Ȃ��`���̓R���p�C���[�Ɉˑ�����̂Œ���
      OPEN (11, FILE=FNAME(1), STATUS='NEW', FORM='UNFORMATTED')
*V�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (12, FILE=FNAME(2), STATUS='NEW', FORM='UNFORMATTED')
*W�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (13, FILE=FNAME(3), STATUS='NEW', FORM='UNFORMATTED')
*P�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (14, FILE=FNAME(4), STATUS='NEW', FORM='UNFORMATTED')
*T�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
      OPEN (15, FILE=FNAME(5), STATUS='NEW', FORM='UNFORMATTED')
*in3d.mac���̃R�����g�s(12-17�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) ITYPE, ICYCLE, NITR, NCYCLE
*in3d.mac���̃R�����g�s(19-20�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) EPSP, OMG
*in3d.mac���̃R�����g�s(22-23�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DT, RE, PR, GR
*in3d.mac���̃R�����g�s(25-26�s��)�̃X�L�b�v
      READ (10,'(A)')
      READ (10,'(A)')
      READ (10,*) DLX, DLY, DLZ, IRELP, METHOD
*�p���̌v�Z�̏ꍇ
      IF (ICYCLE.NE.0) THEN
*U�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (16, FILE=FNAME(6), STATUS='OLD', FORM='UNFORMATTED')
*V�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (17, FILE=FNAME(7), STATUS='OLD', FORM='UNFORMATTED')
*W�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (18, FILE=FNAME(8), STATUS='OLD', FORM='UNFORMATTED')
*P�f�[�^�t�@�C���̃I�[�v��(�����Ȃ��`��)
        OPEN (19, FILE=FNAME(9), STATUS='OLD', FORM='UNFORMATTED')
*T�f�[�^�t�@�C���̃I�[�v��(������ł�T=0.0�̃f�[�^��ǂݍ���)(�����Ȃ��`��)
        OPEN (20, FILE=FNAME(10), STATUS='OLD', FORM='UNFORMATTED')
      END IF
*x�����̊i�q��
      DX  = DLX / FLOAT(NX)
*y�����̊i�q��
      DY  = DLY / FLOAT(NY)
*Z�����̊i�q��
      DZ  = DLZ / FLOAT(NZ)
*�^�����������̊g�U���̌W��(�����ł�Pr)
      VIS = PR
*�G�l���M�[���������̊g�U���̌W��(�����ł�1)
      ALP = 1.0D+0
*���͍��̌W��(�����ł� Gr * Pr**2)
      BUO = GR * PR**2
*������Ȃ畂�͍��̌W���̓[��
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
*       ���͌v�Z�̔����񐔂�NITR�ȏ�̂Ƃ����U�Ƃ݂Ȃ��Čv�Z�I��
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
      CALL TECPLT (FNAME(11))
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      CLOSE (15)
      STOP
      END
*********************************************************************
*                          �����ݒ�
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*�V�K�v�Z�̏ꍇ
      IF ( ICYCLE. EQ. 0 ) THEN
*       U�̏����l�ݒ�
        DO 10 IX = 0,NX
          DO 20 IY = 0,NY+1
            DO 30 IZ = 0,NZ+1
              UN(IX,IY,IZ) = 0.0D0
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
*       V�̏����l�ݒ�
        DO 40 IX = 0,NX+1
          DO 50 IY = 0,NY
            DO 60 IZ = 0,NZ+1
              VN(IX,IY,IZ) = 0.0D0
   60     CONTINUE
   50     CONTINUE
   40   CONTINUE
*       W�̏����l�ݒ�
        DO 70 IX = 0,NX+1
          DO 80 IY = 0,NY+1
            DO 90 IZ = 0,NZ
              WN(IX,IY,IZ) = 0.0D0
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
*       P�̏����l�ݒ�
        DO 100 IX = 0,NX+1
          DO 110 IY = 0,NY+1
            DO 120 IZ = 0,NZ+1
              PO(IX,IY,IZ) = 0.0D0
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE
*-----------------------------------------------------------------------*
* �i���Ӂj���͍��̌v�Z�ŉ��x�̔z����g�p���Ă���̂œ�����ł�T=0�Ƃ��� *
* �������������͐ݒ肷��K�v������D�[���ȊO�̒l������ƕ��͍����v�Z  *
* �����\��������̂Œ��ӁD                                          *
*-----------------------------------------------------------------------*
*       T�̏����l�ݒ�(�̈���͍���(+0.5)�ƒቷ(-0.5)�̒��ԉ��x)
        DO 130 IX = 0,NX+1
          DO 140 IY = 0,NY+1
            DO 150 IZ = 0,NZ+1
              TN(IX,IY,IZ) = 0.0D0
  150     CONTINUE
  140     CONTINUE
  130   CONTINUE
*       T�̋��E�F�E���ǁi��p�jT=-0.5
        DO 160 IY = 0,NY+1
          DO 170 IZ = 0,NZ+1
            TN(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY,IZ)
  170     CONTINUE
  160   CONTINUE
*       T�̋��E�F�����ǁi���M�jT=+0.5
        DO 180 IY = 0,NY+1
          DO 190 IZ = 0,NZ+1
            TN(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - TN(1,IY,IZ)
  190   CONTINUE
  180   CONTINUE
*       T�̋��E�F��ʁi�f�M�j
        DO 200 IX = 1,NX
          DO 210 IZ = 0,NZ+1
            TN(IX,NY+1,IZ) = TN(IX,NY,IZ)
  210     CONTINUE
  200   CONTINUE
*       T�̋��E�F���ʁi�f�M�j
        DO 220 IX = 1,NX
          DO 230 IZ = 0,NZ+1
            TN(IX,0,IZ) = TN(IX,1,IZ)
  230     CONTINUE
  220   CONTINUE
*       T�̋��E�F��ʁi�f�M�j
        DO 240 IX = 1,NX
          DO 250 IY = 1,NY
            TN(IX,IY,NZ+1) = TN(IX,IY,NZ)
  250     CONTINUE
  240   CONTINUE
*       T�̋��E�F�O�ʁi�f�M�j
        DO 260 IX = 1,NX
          DO 270 IY = 1,NY
            TN(IX,IY,0) = TN(IX,IY,1)
  270     CONTINUE
  260   CONTINUE
*�p���v�Z�i���łɂ���v�Z���ʂ���X�^�[�g�j�̏ꍇ
      ELSE
*       U�f�[�^�t�@�C������̓ǂݍ���[Unit No.=16](�����Ȃ��`��)
        READ (16) UN
*       V�f�[�^�t�@�C������̓ǂݍ���[Unit No.=17](�����Ȃ��`��)
        READ (17) VN
*       W�f�[�^�t�@�C������̓ǂݍ���[Unit No.=18](�����Ȃ��`��)
        READ (18) WN
*       P�f�[�^�t�@�C������̓ǂݍ���[Unit No.=19](�����Ȃ��`��)
        READ (19) PO
*---------------------------------------------------------------------*
*    (����) ������̌v�Z�ł�T(=0)�̃t�@�C����ǂݍ��ޕK�v������       *
*---------------------------------------------------------------------*
*       T�f�[�^�t�@�C������̓ǂݍ���[Unit No.=20](�����Ȃ��`��)
        READ (20) TN
        CLOSE (16)
        CLOSE (17)
        CLOSE (18)
        CLOSE (19)
        CLOSE (20)
      END IF
      RETURN
      END
*********************************************************************
*                         ���Ԑi�s
*********************************************************************
      SUBROUTINE ADV
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      TIME   = DT*FLOAT(ICYCLE)
      ICYCLE = ICYCLE + 1
*ICYCLE��100�񖈂ɕ\��
      IF (MOD(ICYCLE,100).EQ.0) THEN
        WRITE (6,2000) ICYCLE
 2000   FORMAT ('  CYC = ',I8)
      END IF
*--------------------------------------------------------------------
* UN -> UO : �K�v�Ȃ����ւ���O��UN��UO����ϓ��ʂ����߂�
* UN : �O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����
* UO : �V�������ԃX�e�b�v�ł̏����l�DUN��ۑ��D
*--------------------------------------------------------------------
      DO 130 IX = 0,NX
        DO 140 IY = 0,NY+1
          DO 150 IZ = 0,NZ+1
            UO(IX,IY,IZ) = UN(IX,IY,IZ)
  150   CONTINUE
  140   CONTINUE
  130 CONTINUE
*--------------------------------------------------------------------
* VN -> VO : �K�v�Ȃ����ւ���O��VN��VO����ϓ��ʂ����߂�
* VN : �O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����
* VO : �V�������ԃX�e�b�v�ł̏����l�DVN��ۑ��D
*--------------------------------------------------------------------
      DO 160 IX = 0,NX+1
        DO 170 IY = 0,NY
          DO 180 IZ = 0,NZ+1
            VO(IX,IY,IZ) = VN(IX,IY,IZ)
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*--------------------------------------------------------------------
* WN -> WO : �K�v�Ȃ����ւ���O��WN��WO����ϓ��ʂ����߂�
* WN : �O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����
* WO : �V�������ԃX�e�b�v�ł̏����l�DWN��ۑ��D
*--------------------------------------------------------------------
      DO 190 IX = 0,NX+1
        DO 200 IY = 0,NY+1
          DO 210 IZ = 0,NZ
            WO(IX,IY,IZ) = WN(IX,IY,IZ)
  210     CONTINUE
  200   CONTINUE
  190 CONTINUE
*--------------------------------------------------------------------
* TN -> TO : �K�v�Ȃ����ւ���O��TN��TO����ϓ��ʂ����߂�
* TN : �O�̎��ԃX�e�b�v�ł̌v�Z�l
* TO : �V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D
*--------------------------------------------------------------------
      DO 220 IX = 0,NX+1
        DO 230 IY = 0,NY+1
          DO 240 IZ = 0,NZ+1
            TO(IX,IY,IZ) = TN(IX,IY,IZ)
  240     CONTINUE
  230   CONTINUE
  220 CONTINUE
      RETURN
      END
*********************************************************************
*                     ���x��̌v�Z
*********************************************************************
      SUBROUTINE CALVEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*--------------------------------------------------------------------
*                       U(IX,IY,IZ)�̌v�Z
*--------------------------------------------------------------------
      DO 10 IX = 1,NX-1
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
*       VV��U(IX,IY,IZ)�ɂ�����V�̕�Ԓl
        VV = (  VO(IX,IY  ,IZ)+VO(IX+1,IY  ,IZ)
     $         +VO(IX,IY-1,IZ)+VO(IX+1,IY-1,IZ) )/4.0D0
*       WW��U(IX,IY,IZ)�ɂ�����W�̕�Ԓl
        WW = (  WO(IX,IY,IZ  )+WO(IX+1,IY,IZ  )
     $         +WO(IX,IY,IZ-1)+WO(IX+1,IY,IZ-1) )/4.0D0
*       �Η���(CNVUX,CNVUY,CNVUZ)���P�����x���㍷���ɂČv�Z
        IF ( UO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVUX = UO(IX,IY,IZ)*( UO(IX,IY,IZ) - UO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVUX = UO(IX,IY,IZ)*( UO(IX+1,IY,IZ) - UO(IX,IY,IZ) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY,IZ) - UO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVUY = VV*( UO(IX,IY+1,IZ) - UO(IX,IY,IZ) ) / DY
        END IF
        IF ( WW. GE. 0.0D0 ) THEN
          CNVUZ = WW*( UO(IX,IY,IZ) - UO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WW. LT. 0.0D0 ) THEN
          CNVUZ = WW*( UO(IX,IY,IZ+1) - UO(IX,IY,IZ) ) / DZ
        END IF
*       x�����̕��͍�(BUOU)�̓[��
        TU = 0.0D0
        BUOU = BUO * TU
*       �g�U��(DIFU)�̌v�Z
        DIFU = VIS*(
     $  +( UO(IX-1,IY,IZ)-2.0D0*UO(IX,IY,IZ)+UO(IX+1,IY,IZ) )/DX**2
     $  +( UO(IX,IY-1,IZ)-2.0D0*UO(IX,IY,IZ)+UO(IX,IY+1,IZ) )/DY**2
     $  +( UO(IX,IY,IZ-1)-2.0D0*UO(IX,IY,IZ)+UO(IX,IY,IZ+1) )/DZ**2
     $              )
*       ���̑��x(U)�̌v�Z
        UN(IX,IY,IZ) = UO(IX,IY,IZ)
     $  + DT*( -CNVUX-CNVUY-CNVUZ+DIFU+BUOU
     $         +(PO(IX,IY,IZ)-PO(IX+1,IY,IZ))/DX )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*--------------------------------------------------------------------
*                    V(IX,IY,IZ)�̌v�Z
*--------------------------------------------------------------------
      DO 40 IX = 1,NX
        DO 50 IY = 1,NY-1
          DO 60 IZ = 1,NZ
*       UU��V(IX,IY,IZ)�ɂ�����U�̕�Ԓl
        UU = (  UO(IX-1,IY  ,IZ)+UO(IX,IY  ,IZ)
     $         +UO(IX-1,IY+1,IZ)+UO(IX,IY+1,IZ) )/4.0D0
*       WW��V(IX,IY,IZ)�ɂ�����W�̕�Ԓl
        WW = (  WO(IX,IY  ,IZ-1)+WO(IX,IY  ,IZ)
     $         +WO(IX,IY+1,IZ-1)+WO(IX,IY+1,IZ) )/4.0D0
*       �Η���(CNVVX,CNVVY,CNVVZ)���P�����x���㍷���ɂČv�Z
        IF ( UU. GE. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX,IY,IZ) - VO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVVX = UU*( VO(IX+1,IY,IZ) - VO(IX,IY,IZ) ) / DX
        END IF
        IF ( VO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVVY = VO(IX,IY,IZ)*( VO(IX,IY,IZ) - VO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVVY = VO(IX,IY,IZ)*( VO(IX,IY+1,IZ) - VO(IX,IY,IZ) ) / DY
        END IF
        IF ( WW. GE. 0.0D0 ) THEN
          CNVVZ = WO(IX,IY,IZ)*( VO(IX,IY,IZ) - VO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WW. LT. 0.0D0 ) THEN
          CNVVZ = WO(IX,IY,IZ)*( VO(IX,IY,IZ+1) - VO(IX,IY,IZ) ) / DZ
        END IF
*       ���͍�(BUOV)�̌v�Z
        TV = ( TO(IX,IY,IZ) + TO(IX,IY+1,IZ) )/2.0D0
        BUOV = BUO*TV
*       �g�U��(DIFV)�̌v�Z
        DIFV = VIS*(
     $  +(VO(IX-1,IY,IZ)-2.0D0*VO(IX,IY,IZ)+VO(IX+1,IY,IZ))/DX**2
     $  +(VO(IX,IY-1,IZ)-2.0D0*VO(IX,IY,IZ)+VO(IX,IY+1,IZ))/DY**2
     $  +(VO(IX,IY,IZ-1)-2.0D0*VO(IX,IY,IZ)+VO(IX,IY,IZ+1))/DZ**2
     $             )
*       ���̑��x(V)�̌v�Z
        VN(IX,IY,IZ) = VO(IX,IY,IZ)
     $  + DT*( -CNVVX-CNVVY-CNVVZ+DIFV+BUOV
     $         +(PO(IX,IY,IZ)-PO(IX,IY+1,IZ))/DY )
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
*--------------------------------------------------------------------
*                     W(IX,IY,IZ)�̌v�Z
*--------------------------------------------------------------------
      DO 70 IX = 1,NX
        DO 80 IY = 1,NY
          DO 90 IZ = 1,NZ-1
*       UU��W(IX,IY,IZ)�ɂ�����U�̕�Ԓl
        UU = (  UO(IX-1,IY,IZ  )+UO(IX,IY,IZ)
     $         +UO(IX-1,IY,IZ+1)+UO(IX,IY,IZ+1) )/4.0D0
*       VV��W(IX,IY,IZ)�ɂ�����V�̕�Ԓl
        VV = (  VO(IX,IY-1,IZ  )+VO(IX,IY,IZ)
     $         +VO(IX,IY-1,IZ+1)+VO(IX,IY,IZ+1) )/4.0D0
*       �Η���(CNVWX,CNVWY,CNVWZ)���P�����x���㍷���ɂČv�Z
        IF ( UU. GE. 0.0D0 ) THEN
          CNVWX = UU*( WO(IX,IY,IZ) - WO(IX-1,IY,IZ) ) / DX
        ELSE IF ( UU. LT. 0.0D0 ) THEN
          CNVWX = UU*( WO(IX+1,IY,IZ) - WO(IX,IY,IZ) ) / DX
        END IF
        IF ( VV. GE. 0.0D0 ) THEN
          CNVWY = VV*( WO(IX,IY,IZ) - WO(IX,IY-1,IZ) ) / DY
        ELSE IF ( VV. LT. 0.0D0 ) THEN
          CNVWY = VV*( WO(IX,IY+1,IZ) - WO(IX,IY,IZ) ) / DY
        END IF
        IF ( WO(IX,IY,IZ). GE. 0.0D0 ) THEN
          CNVWZ = WO(IX,IY,IZ)*( WO(IX,IY,IZ) - WO(IX,IY,IZ-1) ) / DZ
        ELSE IF ( WO(IX,IY,IZ). LT. 0.0D0 ) THEN
          CNVWZ = WO(IX,IY,IZ)*( WO(IX,IY,IZ+1) - WO(IX,IY,IZ) ) / DZ
        END IF
*       ���͍�(BUOW)�̌v�Z
        TW = 0.0D0
        BUOW = BUO * TW
*       �g�U��(DIFW)�̌v�Z
        DIFW = VIS*(
     $  +(WO(IX-1,IY,IZ)-2.0D0*WO(IX,IY,IZ)+WO(IX+1,IY,IZ))/DX**2
     $  +(WO(IX,IY-1,IZ)-2.0D0*WO(IX,IY,IZ)+WO(IX,IY+1,IZ))/DY**2
     $  +(WO(IX,IY,IZ-1)-2.0D0*WO(IX,IY,IZ)+WO(IX,IY,IZ+1))/DZ**2
     $             )
*       ���̑��x(W)�̌v�Z
        WN(IX,IY,IZ) = WO(IX,IY,IZ)
     $  + DT*( -CNVWX-CNVWY-CNVWZ+DIFW+BUOW
     $         +(PO(IX,IY,IZ)-PO(IX,IY,IZ+1))/DZ )
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*���x�̋��E�����̏���
      CALL VELBND
      RETURN
      END
*********************************************************************
*                    ���͏�̌v�Z
*********************************************************************
      SUBROUTINE PRESS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      IXMAX = 0
      IYMAX = 0
      IZMAX = 0
      DMAX = 0.0D0
*P(IX,IY,IZ)�̌v�Z
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
          DEL = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          DIV = ( UN(IX,IY,IZ) - UN(IX-1,IY  ,IZ  ) )/DX
     $        + ( VN(IX,IY,IZ) - VN(IX  ,IY-1,IZ  ) )/DY
     $        + ( WN(IX,IY,IZ) - WN(IX  ,IY  ,IZ-1) )/DZ
          IF ( DABS(DIV). GE. DABS(DMAX) ) THEN
            IXMAX = IX
            IYMAX = IY
            IZMAX = IZ
            DMAX = DIV
          END IF
          DELP = - OMG * DIV / DEL
          PO(IX,  IY  ,IZ  ) = PO(IX,  IY,  IZ  ) + DELP
          UN(IX,  IY  ,IZ  ) = UN(IX,  IY,  IZ  ) + DT/DX*DELP
          UN(IX-1,IY  ,IZ  ) = UN(IX-1,IY,  IZ  ) - DT/DX*DELP
          VN(IX,  IY  ,IZ  ) = VN(IX,  IY,  IZ  ) + DT/DY*DELP
          VN(IX,  IY-1,IZ  ) = VN(IX,  IY-1,IZ  ) - DT/DY*DELP
          WN(IX,  IY  ,IZ  ) = WN(IX,  IY,  IZ  ) + DT/DZ*DELP
          WN(IX,  IY  ,IZ-1) = WN(IX,  IY,  IZ-1) - DT/DZ*DELP
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*
* ���͂̑��ΐ��Ɋւ��鏈��(IRELP=1�Ȃ�ȉ��̏������s��)
      IF (IRELP.EQ.1) THEN
        POSTN = PO(1,1,1)
        DO 40 IX = 1,NX
          DO 50 IY = 1,NY
            DO 60 IZ = 1,NZ
              PO(IX,IY,IZ) = PO(IX,IY,IZ)-POSTN
   60       CONTINUE
   50     CONTINUE
   40   CONTINUE
      END IF
*
* IFLG=1�Ȃ�A���̎��𖞂����Ă��Ȃ��Ɣ��肵�Ăш��͌v�Z���s���D
      IF ( DABS(DMAX). GE. EPSP ) IFLG = 1
* ���͌v�Z�̉񐔂�100�񂲂Ƃɕ\��
      IF (MOD(ITR,100).EQ.0) THEN
        WRITE (6,2000) ITR, IXMAX, IYMAX, IZMAX, DMAX
 2000   FORMAT (' Iteration=',I8,'   Div(max)(',3I6,')=',1PE13.5)
      END IF
*�V���ɓ���ꂽ���x��p���ċ��E��������������
      CALL VELBND
      RETURN
      END
*********************************************************************
*                     ���x��̌v�Z
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*T(IX,IY,IZ)�̌v�Z
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
*         UUT,VVT,WWT�͂��ꂼ��T(IX,IY,IZ)�ɂ�����U,V�̕�Ԓl
          UUT = ( UO(IX,IY,IZ) + UO(IX-1,IY,  IZ  ) ) / 2.0D0
          VVT = ( VO(IX,IY,IZ) + VO(IX  ,IY-1,IZ  ) ) / 2.0D0
          WWT = ( WO(IX,IY,IZ) + WO(IX  ,IY,  IZ-1) ) / 2.0D0
*         �Η���(CNVTX,CNVTY,CNVTZ)���P�����x���㍷���ɂČv�Z
          IF ( UUT. GE. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX,IY,IZ) - TO(IX-1,IY,IZ) ) / DX
          ELSE IF ( UUT. LT. 0.0D0 ) THEN
            CNVTX = UUT*( TO(IX+1,IY,IZ) - TO(IX,IY,IZ) ) / DX
          END IF
          IF ( VVT. GE. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY,IZ) - TO(IX,IY-1,IZ) ) / DY
          ELSE IF ( VVT. LT. 0.0D0 ) THEN
            CNVTY = VVT*( TO(IX,IY+1,IZ) - TO(IX,IY,IZ) ) / DY
          END IF
          IF ( WWT. GE. 0.0D0 ) THEN
            CNVTZ = WWT*( TO(IX,IY,IZ) - TO(IX,IY,IZ-1) ) / DZ
          ELSE IF ( WWT. LT. 0.0D0 ) THEN
            CNVTZ = WWT*( TO(IX,IY,IZ+1) - TO(IX,IY,IZ) ) / DZ
          END IF
*         �g�U��(DIFT)�̌v�Z
          DIFT = ALP*(
     $    +( TO(IX-1,IY,IZ)-2.0D0*TO(IX,IY,IZ)+TO(IX+1,IY,IZ) )/DX**2
     $    +( TO(IX,IY-1,IZ)-2.0D0*TO(IX,IY,IZ)+TO(IX,IY+1,IZ) )/DY**2
     $    +( TO(IX,IY,IZ-1)-2.0D0*TO(IX,IY,IZ)+TO(IX,IY,IZ+1) )/DZ**2
     $                )
*         ���̎��Ԃ�T�̌v�Z
          TN(IX,IY,IZ) = TO(IX,IY,IZ) + DT*( -CNVTX-CNVTY-CNVTZ+DIFT )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*���E�����̏���
      CALL TBND
      RETURN
      END
*********************************************************************
*                 ���x�̋��E�����̏���
*********************************************************************
      SUBROUTINE VELBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*U�i�E�ʁj
      DO 10 IY = 1,NY
        DO 20 IZ = 1,NZ
          UN(NX,IY,IZ) = 0.0D0
   20   CONTINUE
   10 CONTINUE
*U�i���ʁj
      DO 30 IY = 1,NY
        DO 40 IZ = 1,NZ
          UN(0,IY,IZ) = 0.0D0
   40   CONTINUE
   30 CONTINUE
*U�i��ʁj
      DO 50 IX = 0,NX
        DO 60 IY = 1,NY
          UN(IX,IY,NZ+1) = -UN(IX,IY,NZ)
   60   CONTINUE
   50 CONTINUE
*U�i�O�ʁj
      DO 70 IX = 0,NX
        DO 80 IY = 1,NY
          UN(IX,IY,0) = -UN(IX,IY,1)
   80   CONTINUE
   70 CONTINUE
*U�i��ʁj
      DO 90 IX = 0,NX
        DO 100 IZ = 0,NZ+1
          UN(IX,NY+1,IZ) = -UN(IX,NY,IZ)
  100   CONTINUE
   90 CONTINUE
*U�i���ʁj
      DO 110 IX = 0,NX
        DO 120 IZ = 0,NZ+1
          UN(IX,0,IZ) = -UN(IX,1,IZ)
  120   CONTINUE
  110 CONTINUE
*V�i�E�ʁj
      DO 200 IY = 1,NY-1
        DO 210 IZ = 1,NZ
          VN(NX+1,IY,IZ) = -VN(NX,IY,IZ)
  210   CONTINUE
  200 CONTINUE
*V�i���ʁj
      DO 220 IY = 1,NY-1
        DO 230 IZ = 1,NZ
          VN(0,IY,IZ) = -VN(1,IY,IZ)
  230   CONTINUE
  220 CONTINUE
*V�i��ʁj
      DO 240 IX = 0,NX+1
        DO 250 IZ = 0,NZ+1
          VN(IX,NY,IZ) = 0.0D0
  250   CONTINUE
  240 CONTINUE
*V�i���ʁj
      DO 260 IX = 0,NX+1
        DO 270 IZ = 0,NZ+1
          VN(IX,0,IZ) = 0.0D0
  270   CONTINUE
  260 CONTINUE
*V�i��ʁj
      DO 280 IX = 0,NX+1
        DO 290 IY = 1,NY-1
          VN(IX,IY,NZ+1) = -VN(IX,IY,NZ)
  290   CONTINUE
  280 CONTINUE
*V�i�O�ʁj
      DO 300 IX = 0,NX+1
        DO 310 IY = 1,NY-1
          VN(IX,IY,0) = -VN(IX,IY,1)
  310   CONTINUE
  300 CONTINUE
*W�i�E�ʁj
      DO 410 IY = 1,NY
        DO 420 IZ = 1,NZ-1
          WN(NX+1,IY,IZ) = -WN(NX,IY,IZ)
  420   CONTINUE
  410 CONTINUE
*W�i���ʁj
      DO 430 IY = 1,NY
        DO 440 IZ = 1,NZ-1
          WN(0,IY,IZ) = -WN(1,IY,IZ)
  440   CONTINUE
  430 CONTINUE
*W�i��ʁj
      DO 450 IX = 0,NX
        DO 460 IY = 0,NY+1
          WN(IX,IY,NZ) = 0.0D0
  460   CONTINUE
  450 CONTINUE
*W�i�O�ʁj
      DO 470 IX = 0,NX
        DO 480 IY = 0,NY+1
          WN(IX,IY,0) = 0.0D0
  480   CONTINUE
  470 CONTINUE
*W�i��ʁj
      DO 490 IX = 0,NX+1
        DO 500 IZ = 1,NZ-1
          WN(IX,NY+1,IZ) = -WN(IX,NY,IZ)
  500   CONTINUE
  490 CONTINUE
*W�i���ʁj
      DO 510 IX = 0,NX+1
        DO 520 IZ = 1,NZ-1
          WN(IX,0,IZ) = -WN(IX,1,IZ)
  520   CONTINUE
  510 CONTINUE
      RETURN
      END
*********************************************************************
*                 ���x�̋��E�����̏���
*********************************************************************
      SUBROUTINE TBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
*�E��
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          TN(NX+1,IY,IZ) = 2.0D0 * ( -0.5D0 ) - TN(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*����
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          TN(0,IY,IZ) = 2.0D0 * ( +0.5D0 ) - TN(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*���
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          TN(IX,NY+1,IZ) = TN(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*����
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          TN(IX,0,IZ) = TN(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*���
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          TN(IX,IY,NZ+1) = TN(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*�O��
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          TN(IX,IY,0) = TN(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                        �f�[�^�o��
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
       WRITE (11) UN
       WRITE (12) VN
       WRITE (13) WN
       WRITE (14) PO
       WRITE (15) TN
      RETURN
      END
*********************************************************************
*                  Tecplot�p�f�[�^�o��
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D3 / VIS,ALP,BUO
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
      COMMON / D6 / DMAX
      COMMON / D7 / ITYPE
      COMMON / ARRAY1 /
     $ UO(0:NX0,  0:NY0+1,0:NZ0+1),UN(0:NX0  ,0:NY0+1,0:NZ0+1),
     $ VO(0:NX0+1,0:NY0  ,0:NZ0+1),VN(0:NX0+1,0:NY0  ,0:NZ0+1),
     $ WO(0:NX0+1,0:NY0+1,0:NZ0  ),WN(0:NX0+1,0:NY0+1,0:NZ0  ),
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
       CHARACTER*20 FNAME
*
      OPEN (21,FILE=FNAME,STATUS='NEW')
      WRITE (21,*) 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "T"'
      NX1 = NX+1
      NY1 = NY+1
      NZ1 = NZ+1
      WRITE (21,4000) NX1,NY1,NZ1
 4000 FORMAT (1H ,'ZONE I=',I3,',J=',I3,',K=',I3,', F=POINT')
      DO 10 IZ = 0,NZ
         DO 20 IY = 0,NY
           DO 30 IX = 0,NX
             X = DX * FLOAT(IX)
             Y = DY * FLOAT(IY)
             Z = DZ * FLOAT(IZ)
             U = ( UN(IX  ,IY,  IZ  )+UN(IX  ,IY+1,IZ  )
     $            +UN(IX  ,IY,  IZ+1)+UN(IX  ,IY+1,IZ+1) )/4.0D0
             V = ( VN(IX  ,IY,  IZ  )+VN(IX+1,IY  ,IZ  )
     $            +VN(IX  ,IY,  IZ+1)+VN(IX+1,IY,  IZ+1) )/4.0D0
             W = ( WN(IX  ,IY  ,IZ  )+WN(IX  ,IY+1,IZ  )
     $            +WN(IX+1,IY  ,IZ  )+WN(IX+1,IY+1,IZ  ) )/4.0D0
             T = ( TN(IX  ,IY  ,IZ  )+TN(IX+1,IY  ,IZ  )
     $            +TN(IX  ,IY+1,IZ  )+TN(IX+1,IY+1,IZ  )
     $            +TN(IX  ,IY  ,IZ+1)+TN(IX+1,IY  ,IZ+1)
     $            +TN(IX  ,IY+1,IZ+1)+TN(IX+1,IY+1,IZ+1) )/8.0D0
            WRITE (21,4010) X,Y,Z,U,V,W,T
 4010       FORMAT (1H ,7(1PE11.3))
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      CLOSE(21)
      RETURN
      END
