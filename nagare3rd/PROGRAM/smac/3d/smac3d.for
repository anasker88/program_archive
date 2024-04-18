*********************************************************************
* �t�@�C����  �FSMAC3D.FOR                                          *
* �^�C�g��    �FSMAC�@�ɂ��3�����M������̓v���O����               *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� �o�C�I�E���p���w��              *
* �����      �F2011.11.01                                          *
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
*    ITYPE......1:�����v�Z 2:�񓙉��v�Z(�P�R�C�P�S�s���Q��)         *
*    ICYCLE.....�v�Z�J�n�̃T�C�N�����i����T=ICYCLE*DT�j             *
*            0�Ȃ�V�K�Ńv���O�����ɂ��鏉�������ɂ��������Čv�Z�J�n*
*            0�ȊO�̒l�Ȃ�p���v�Z                                  *
*            ���� TIME = ICYCLE * DT                                *
*    NITR....���͕␳�̐��`�V�X�e����@(�����@�E�N�����t������Ԗ@) *
*            �̂��߂̍ő唽����                                   *
*    NCYCLE..�v�Z�I���T�C�N����                                     *
*    EPSP....��������l                                             *
*            ���͕␳�̐��`�V�X�e����@(�����@,�N�����t������Ԗ@)  *
*            �̎����]���Ɏg�p                                       *
*    OMG....���͕␳�̂��߂̊ɘa�W��: (point, line-)SOR�@�̉����W�� *
*    DT.....���ԍ���                                                *
*    RE.....���C�m���Y��, PR.....�v�����g����, GR.....�O���X�z�t��  *
*    DLX.....��͗̈�̉���(DLX=NX*DX):DX�͈��œ��Ԋu�i�q         *
*    DLY.....��͗̈�̍���(DLY=NY*DY):DY�͈��œ��Ԋu�i�q         *
*            �i�q���i���Ԋu�jDX,DY�̓v���O�����̒��ŋ��߂�D        *
*   IRELP...���͂̊�l�Ɖ���1���Ɨ����̐ݒ�                       *
*           0 : �s��Ȃ�                                            *
*             <SMAC�ɂ����鈳�͕␳�̐��`�V�X�e����@�Ɋւ���>      *
*             -> 1���]���ȉ���1�����߂�̂�                       *
*             -> ���ږ@�ł͓��ٍs��̖��ɑ���. �w�I�ɂ͐������Ȃ�.*
*                �ۂߌ덷���Ȃ���Ή��͓����Ȃ��D                 *
*           1 : ���͊�𔽉f����1���Ɨ��ȉ������߂�D             *
*           2 : ���͊��ݒ� P(1,1)=0                             *
*               1���]���̉���1������,���͊�lP(1,1)=0��ݒ肷��.*
*         (IRELP=0�̌v�Z��PD�����߁CPD(1,1)����������P(1,1)=0�Ƃ���)*
*    METHOD..���͕␳�̐��`�V�X�e����@�̃A���S���Y��               *
*            ���ׂăo���h�}�g���b�N�X�p�ɍœK�����Ă���             *
*            1: ���ږ@ -> �K�E�X�̏����@                            *
*               IRELP=1�Ƃ���K�v������D(�ۂߌ덷�ɂ��CIREP=0,2  *
*               �Ƃ��Ă����𓾂���ꍇ�����邪���w�I�ɐ������Ȃ�.)*
*            2: �����@1 -> point-SOR �@                             *
*            3: �����@2 -> line-SOR �@                              *
*               OMG=1�ŏ\���D�傫����������Ɣ��U����               *
*            4: �N�����t������Ԗ@1 -> �����c���@                   *
*            5: �N�����t������Ԗ@2 -> Bi-CGSTAB�@                  *
*    �p�����[�^�t�@�C���̐��l�ɂ��āC                             *
*    FORTRAN �v���O����.....�ł���Δ{���x�����ŗ^����"1.0D0,1.0d0" *
*            C�v���O�����Ƌ��p������"1.0E0,1.0e0"�Ƃ��Ă����͂Ȃ� *
*    C �v���O����....."1.0E0,1.0e0"�Ƃ��ė^����                     *
*    TECPLOT�p�f�[�^�����������o�̓t�@�C���͏����Ȃ��`���ŁC�g�p    *
*    ����R���p�C���[�Ɉˑ�����D�R���p�C���[�Ɉˑ����Ȃ��`���ɂ��� *
*    �ɂ́C�e�ʂ͑����邪�����t���`���ɂ���΂悢�D                 *
*                                                                   *
*  �����͂̑��ΐ��ɂ��āi�{���F14.5.3�����Q�Ɓj                   *
*    ���͕␳�̐��`�V�X�e���������ɂ�����C�{���̂悤�ȁC���E�ɂ� *
*    ���鑬�x�����m�̏ꍇ�C�W���s�񂪓���(matrix singular)�ƂȂ�C  *
*    ��������͈ꎟ�]���ƂȂ�C�����̉������݂��邱�ƂƂȂ�D     *
*    ����ɑ΂������Ƃ��āC                                       *
*    (1) ���ږ@��p����Ƃ��́C�������̐������炷�K�v������D       *
*    �T�u���[�`��PRESS�ɂ����āC��Ƃ��� ��P(1,1)=P(1,1)=0 �ƌŒ� *
*    �ł���悤�ɂȂ��Ă���D                                       *
*    (2) �N�����t������Ԗ@���܂ޔ����@�ɂ����ẮC�ӎ����Ȃ��Ă悢 *
*    �ꍇ�������D                                                   *
*    ���w�I�Ɍ����ɂ����΁C�ꎟ�]���ȉ��̂�����1��������Ƃ���  *
*    ���Ƃ́C�������Ƃ͂����Ȃ����C���͂̒l�́C��Βl�ł͂Ȃ����Βl *
*    �݂̂����ƂȂ�D                                             *
*    �������C�v�Z��i�s�����Ă䂭�ɂ�C���͂̐�Βl���傫���Ȃ��� *
*    �䂭�悤�ȏꍇ�́C�ꎟ�]���ȉ��̂�����1�����߂Ă���C��l��*
*    �ݒ肵���ق����悢�D�͂��߂���C��l��ݒ肵�Ĉꎟ�Ɨ��ȉ��� *
*    ����ɂ́C�ʏ�C�v�Z���Ԃ������Ȃ�D                           *
*                                                                   *
*  ���ϐ��E�z��̐���                                               *
*    ICYCLE  -----> ���Ԑi�s�̂��߂̃J�E���^                        *
*    ITR     -----> ���͕␳�v�Z�̂��߂̔����񐔂̃J�E���^          *
*                 �i�N�����t������Ԗ@���܂ޔ����@�ɂ����Ďg�p�j    *
*    IX,IY,IZ   -----> ��̐}���Q��                                 *
*    UO,VO,WO,TO-----> ���͕␳�v�Z���s���O�̒l                     *
*    UN,VN,WN,TN-----> �V���Ȉ��͕␳��p���Čv�Z���ꂽ�l           *
*    PD         -----> ���͕␳                                     *
*                                                                   *
*  ��SMAC�@�ɂ�������`�V�X�e����@�ɂ���                         *
*    �Î~��Ԓ���ȂǁC���x�ꂪ�[���ł������肷��ƁCMETHOD=1,4,5�� *
*    �����ĉ��𓾂��Ȃ��ꍇ������D                               *
*    METHOD=1: �ꎟ�Ɨ��ȉ������߂�̂����w�I�ɐ������̂�IRELP=1��  *
*              �Ƃ��邱�ƁD�R���p�C���[�ɂ���ẮCIRELP=0,2�Ƃ���  *
*              ���ۂߌ덷�Ȃǂɂ��C����������ꍇ������D       *
*    METHOD=4: �����l�x�N�g��\vec{X}_{0}���[���ƂȂ����肷��ƁC�T��*
*              �x�N�g���̌v�Z���[�����Z�ƂȂ�̂œK���ɏ����l��ݒ� *
*              ����K�v������D�{�v���O�����ł́C�����l�Ƃ��čŐV�� *
*              �l��p����悤�ɂ��āC�\�Ȍ���J��Ԃ��񐔂����Ȃ� *
*              �Ȃ邱�Ƃ�D��Ƃ��Ă���C�[�����Z�̂Ƃ��� POINT-SOR *
*              �@�ɐ؂�ւ���悤�ɂ��Ă���̂ŁC���̂Ƃ��́COMG��  *
*              �l��p���邱�ƂƂȂ�D                               *
*    METHOD=5: METHOD=4�ɓ����D                                     *
*                                                                   *
*    �{�v���O�����́CHAMAC�@�̃v���O����hsmac3d.for�����Ƃɂ��č쐬 *
*    ���Ă���DHSMAC��SMAC���{���I�ɓ����ł���̂ŁC�v���O�����̕ύX*
*    �́C���͕␳���j���[�g���@�ōs��(HSMAC�@)����ɁC���`�V�X�e��*
*    ��@�ōs���悤�ɂ���΂悢�D�Ȃ��C�ł������Chsmac3d.for�Ƃ̈�*
*    �����킩��₷�����邽�߁CSMAC�@�ɂ����ĐV���ɕt������������ *
*    �́C"--- SMAC ---"�Ƃ��ăR�����g��}�����Ă���D               *
*                                                                   *
*********************************************************************
      PROGRAM SMAC3D
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     ���͕␳�̐��`�V�X�e���p
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),A6(NE0),
     $          A7(NE0),B(NE0),X(NE0)
*     ��Ɨp�z��
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),W6(NE0),
     $          W7(NE0),W8(NE0),W9(NE0)
*     �o���h�K�E�X�����@�̂��߂̌W���s��z��
      DIMENSION AGB(-NXY0:NXY0,NE0)
*
      CHARACTER FNAME(12)*20
*
*x�����̊i�q������
      NX  = NX0
*y�����̊i�q������
      NY  = NY0
*Z�����̊i�q������
      NZ  = NZ0
*�������̐�(���͕␳�Ɋւ��関�m���̐�) --- SMAC ---
      NE  = NE0
* NX * NY
      NXY=NXY0
*�p�����[�^�t�@�C���̃I�[�v��
      OPEN (10,FILE='in3d.mac',STATUS='OLD')
*�o�̓t�@�C�����̓ǂݍ���
      DO 10 I = 1,11
        READ (10,'(A20)') FNAME(I)
   10 CONTINUE
*U�̌v�Z���ʏo�͗p�t�@�C���I�[�v��(�����Ȃ��`��)
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
      READ (10,*) ITYPE, ICYCLE0, NITR, NCYCLE
      ICYCLE = ICYCLE0
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
*���͍��̌W��( �����ł� Pr * Ra = Pr * ( Pr * Gr) )
      BUO = ( GR * PR**2 )
*������Ȃ畂�͍��̌W���̓[���ɐݒ�
      IF (ITYPE.EQ.1) BUO = 0.0D0
*�����l�̐ݒ�
      CALL CINITI
      DO 20 I=1,NE
        A1(I)=0.0D0
        A2(I)=0.0D0
        A3(I)=0.0D0
        A4(I)=0.0D0
        A5(I)=0.0D0
        A6(I)=0.0D0
        A7(I)=0.0D0
        B(I)=0.0D0
        X(I)=0.0D0
        W1(I)=0.0D0
        W2(I)=0.0D0
        W3(I)=0.0D0
        W4(I)=0.0D0
        W5(I)=0.0D0
        W6(I)=0.0D0
        W7(I)=0.0D0
        W8(I)=0.0D0
        W9(I)=0.0D0
        DO 25 II = -NXY,NXY
          AGB(II,I)=0.0D0
   25   CONTINUE
   20 CONTINUE
*���Ԑi�s�̂��߂̖߂�_
  700 CONTINUE
*���Ԑi�s
      CALL ADV
*���͕␳�̌v�Z�������������ǂ����̃p�����[�^(IFLG)��������(�����@�ł̂ݗL��)
*IFLG -> 0:���� 1:���U(�ݒ肳�ꂽ���e��NITR�ȉ��ŉ��������Ȃ�)
      IFLG = 0
C     --- SMAC ---
C     �{�v�Z�ɂ����ẮCICYCLE=0�̂Ƃ��̏��������͑��x��ƈ��͏���[����
C     ���Ă���̂ŁC���͕␳�̐��`�V�X�e���������̂�����ƂȂ邽�߁C
C     ICYCLE=1�̂Ƃ������͉��x�̌v�Z�փW�����v����悤�ɂ���D
C     �v�Z��������ɉ����ēK�X�폜���邢�͕ύX
      IF (ICYCLE.EQ.1) GOTO 720
*���x��̌v�Z
      CALL CALVEL
*���͏�̌v�Z
C     CALL PRESS
      CALL PRESS (A1,A2,A3,A4,A5,A6,A7,B,X,
     $            AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9,NXY)
C     ���x��ƈ��͏ꂪ�[���ƂȂ�Ƃ��̃X�L�b�v��--- SMAC ---
  720 CONTINUE
*     ���͏�̌v�Z�����������Ƃ�
      IF ( IFLG. EQ. 0 ) THEN
*       �񓙉���v�Z�̏ꍇ
        IF ( ITYPE.EQ.2 ) THEN
*         ���x����v�Z
          CALL CALTEM
        END IF
*     ���͏�̌v�Z���������Ă��Ȃ��Ƃ�
      ELSE IF ( IFLG. EQ. 1 ) THEN
        WRITE (6,*) ' NOT CONVERGE ! '
C       �f�[�^���o�͂��ċ����I��
        CALL PROUT
        GO TO 900
      END IF
*     ���Ԑi�s�J�E���^(ICYCLE)��NCYCLE��菬������
      IF ( ICYCLE. LT. NCYCLE ) THEN
        GO TO 700
*     ���Ԑi�s�T�C�N����NCYCLE�ɂȂ�����v�Z�I��
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
*
*********************************************************************
*                          �����ݒ�
*********************************************************************
      SUBROUTINE CINITI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
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
C           --- SMAC ---
            PD(IX,IY,IZ) = 0.0D0
            PN(IX,IY,IZ) = 0.0D0
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
C     --- SMAC ---
        READ (19) PN,PD
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
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
* VO : �V�������ԃX�e�b�v�ł̏����l�CVN��ۑ��D
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
* WO : �V�������ԃX�e�b�v�ł̏����l�CVN��ۑ��D
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
* TN : �O�̎��ԃX�e�b�v�ł̒l
* TO : �V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D
* PN -> PO : �K�v�Ȃ����ւ���O��PN��PO����ϓ��ʂ����߂�
* PN : �O�̎��ԃX�e�b�v�ł̒l
* PO : �V�������ԃX�e�b�v�ł̏����l�DPN��ۑ��D
*--------------------------------------------------------------------
      DO 220 IX = 0,NX+1
        DO 230 IY = 0,NY+1
          DO 240 IZ = 0,NZ+1
            TO(IX,IY,IZ) = TN(IX,IY,IZ)
C           --- SMAC ---
            PO(IX,IY,IZ) = PN(IX,IY,IZ)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
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
*                                                                   *
*�m���͕␳�̐��`�V�X�e���Ɋւ���z��̐����n                       *
*                                                                   *
* ���͕␳(PD)�̐��`�V�X�e���� A_{i,j} PD_{i} = B_{i} �Ƃ���D      *
* �L�������ߎ���p���ė��U�����Ă��邱�Ƃ��疾�炩�Ȃ悤�ɁC        *
* 1. A_{i,j}�̑啔���̓[��                                          *
* 2. ��[���v�f�͑a�ł���C�K���I�ɕ���ł���                       *
* ���������āCA_{i,j}���ׂĂ��L��������̂͌v�Z�e�ʂ̃��_�ł���̂� *
* �{�v���O�����ł́C�ȉ��̂悤�ȋK���ɂ��������ċL������D          *
* 3�������ł���΁C���͕␳�Ɋւ���|�A�\���������ɂ����āC2�K��  *
* �������͋ߗׂ�7�_�ŕ\����D                                       *
*                                                                   *
* ���p����z��(1�����z��Ɋi�[)�F�N�����t������Ԗ@���܂ޔ����@     *
*   �W���s��     ---> A1, A2, A3, A4, A5, A6, A7                    *
*   ���m�x�N�g�� ---> B                                             *
*   ���m�x�N�g�� ---> X                                             *
*   ���m���̐�(���m��PD�̐�)   ---> NX * NY * NZ                    *
*                                                                   *
*    A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)               *
*                                                                   *
*    i=3, j=3, k=3 (m=59): PD(3,3,3) �ɂ����� A1,A2,A3,A4,A5,A6,A7  *
*                                                                   *
*        +------------------/-----+                                 *
* (NY)5  |    |    |    |  / |    |                                 *
*        +---------------A7-------+    �ʂ��ԍ��̋K��               *
*     4  |    |    |A4  /    |    |    m=(k-1)*(NX*NY)+(i-1)*NY+j   *
*        +-------------/----------+                                 *
*     3  |    | A1 |A3  | A5 |    |    A1(m)                        *
*        +---------/--------------+    A2(m)                        *
*     2  |    |   /|A2  |    |    |    A3(m)                        *
*        +------A6----------------+    A4(m)                        *
*     1  |      /  |    |    |    |    A5(m)                        *
*        +-----/------------------+    A6(m)                        *
*     j   i=1 /   2    3    4    5(NX) A7(m)                        *
*            /                         B(m)                         *
*         k=1,5(NZ)                    NXY = NX * NY                *
*                                                                   *
*    B(NE) : ���m�x�N�g��                                           *
*                                                                   *
*    X(NE) : ���m�x�N�g��(�e�T�u���[�`���ł��ꂪ���܂�)             *
*                                                                   *
* ���p����z��F�K�E�X�̏����@�ɂ�钼�ږ@(�o���h�s��Ƃ��Ċi�[)    *
*    NXY = NX * NY                                                  *
*    A(-NXY:NXY,NE) : 3���������ߎ��ɂ��K���I��Ώ̍s��           *
*          -NXY:NXY->��̐}�ōl�����M=59�̂Ƃ�A1����A7�܂ł�m��    *
*                   A1��m�� M"-NY"                                  *
*                   A2��m�� M"-1 "                                  *
*                   A3��m�� M"+0 "                                  *
*                   A4��m�� M"+1 "                                  *
*                   A5��m�� M"+NY"                                  *
*                   A6��m�� M"-(NX*NY)"                             *
*                   A7��m�� M"+(NX*NY)"                             *
*              ���̂悤�ɁCm���Œ肵���Ƃ��C" "�ň͂܂ꂽ7�̒l    *
*              (-NY,-1,0,1,NY,-NXY,NXY)                             *
*     NE -> ��q�̂悤�� m �͑S���� NE ��`�����                 *
*                                                                   *
* �����`�V�X�e�����\������ۂ̒���                                  *
* ���`�V�X�e�����\������ہC���E�������ǂ��Ŕ��f�����邩�ɂ���āC  *
* �v���O���~���O���قȂ�D�{�v���O�����ł́C���E�����́C�W���s���  *
* �쐬����ۂɔ��f�����C���`�V�X�e���̉�@�ɂ����Ă͐�� AX=B �̂�  *
* �ɒ��ڂ��闧����Ƃ�D                                            *
* ��q�̂悤�ȗ���Ƃ͈قȂ�C���`�V�X�e���̉�@�ɂ����ċ��E������  *
* ���f�����邱�Ƃ��ł��邪�C�����ł́C���`�V�X�e���̉�@�̃T�u���[  *
* �`���ɔėp�����������邱�Ƃ�D�悳�����D                          *
*                                                                   *
*********************************************************************
      SUBROUTINE PRESS (A1,A2,A3,A4,A5,A6,A7,B,X,
     $                  AGB,W1,W2,W3,W4,W5,W6,W7,W8,W9,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     ���͕␳�̐��`�V�X�e���p
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
*     ��Ɨp�z��
      DIMENSION W1(NE0),W2(NE0),W3(NE0),W4(NE0),W5(NE0),
     $          W6(NE0),W7(NE0),W8(NE0),W9(NE0)
*     �o���h�K�E�X�����@�̂��߂̌W���s��z��
      DIMENSION AGB(-NXY:NXY,NE0)
*
*���`�V�X�e���֌W�̔z��̏�����
      DO 1 IX = 1,NX
        DO 2 IY = 1,NY
          DO 3 IZ = 1,NZ
          K = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          B(K) = 0.0D0
          A1(K) = 0.0D0
          A2(K) = 0.0D0
          A3(K) = 0.0D0
          A4(K) = 0.0D0
          A5(K) = 0.0D0
          A6(K) = 0.0D0
          A7(K) = 0.0D0
          X(K)  = 0.0D0
    3     CONTINUE
    2   CONTINUE
    1 CONTINUE
*P(IX,IY,IZ)�̌v�Z
      DO 10 IX = 1,NX
        DO 20 IY = 1,NY
          DO 30 IZ = 1,NZ
          DIV(IX,IY,IZ) = ( UN(IX,IY,IZ) - UN(IX-1,IY  ,IZ  ) )/DX
     $                  + ( VN(IX,IY,IZ) - VN(IX  ,IY-1,IZ  ) )/DY
     $                  + ( WN(IX,IY,IZ) - WN(IX  ,IY  ,IZ-1) )/DZ
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
*     �W���s��̍쐬
*         �v�Z�̈��1�_����ɓ����̓_�i���z�Z�����܂߂čl�����2�_�����������̓_�j
*         �Ɋւ��ČW���s����쐬�D�c��͋��E�����𔽉f������ PDBNDC �Őݒ肷��
      DO 40 IX = 2,NX-1
        DO 50 IY = 2,NY-1
          DO 60 IZ = 2,NZ-1
            I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
            A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
            A4(I) = - 1.0D0 / DY**2 * DT
            A2(I) = - 1.0D0 / DY**2 * DT
            A1(I) = - 1.0D0 / DX**2 * DT
            A5(I) = - 1.0D0 / DX**2 * DT
            A6(I) = - 1.0D0 / DZ**2 * DT
            A7(I) = - 1.0D0 / DZ**2 * DT
            B(I) = -DIV(IX,IY,IZ)
            X(I) = PD(IX,IY,IZ)
   60     CONTINUE
   50   CONTINUE
   40 CONTINUE
* --- SMAC --- ���E�������W���s��ɔ��f������
      CALL PDBNDC (A1,A2,A3,A4,A5,A6,A7,B)
* --- SMAC --- ���͕␳ P'(PD) �Ɋւ���|�A�\���������̉�@
*�T�u���[�`���Ɉ�ʐ����������邽�߁CNX,NY,NZ,NE�������Ƃ��ēn��
*COMMON���Œ�`����Ă���l�ł�����̂ŁC���̂܂ܓn���Ȃ��I
      NNX = NX
      NNY = NY
      NNZ = NZ
      NNE = NE
      NXY=NX*NY
C   1. ���ږ@ : �o���h�}�g���b�N�X�ɂ��K�E�X�̏����@
      IF (METHOD.EQ.1)
     $  CALL GB      (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                AGB)
C   2. �����@1 : point-SOR �@
      IF (METHOD.EQ.2)
     $  CALL PSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY)
C   3. �����@2 : line-SOR �@ : OMG=1�ŏ\���D�傫����������Ɣ��U����
      IF (METHOD.EQ.3)
     $  CALL LSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6,W7,W8,W9)
C   4. �N�����t������Ԗ@1 : �����c���@
      IF (METHOD.EQ.4)
     $  CALL CRB     (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6)
C   5. �N�����t������Ԗ@2 : Bi-CGSTAB
      IF (METHOD.EQ.5)
     $  CALL BICGB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY,
     $                W1,W2,W3,W4,W5,W6,W7)
C   METHOD=4,5�̂Ƃ��̒T���x�N�g���v�Z���̃[�����Z�΍�
      IF (IFLG.EQ.2)
     $  CALL PSORB   (A1,A2,A3,A4,A5,A6,A7,B,X,NNX,NNY,NNZ,NNE,NXY)
      DO 70 IX = 1,NX
        DO 80 IY = 1,NY
          DO 90 IZ = 1,NZ
          K = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
*         ���͂̑��ΐ��̏��� : ��l�̐ݒ�
*         1���Ɨ��ȉ������߂�ꍇ�͈��͂̊�_�������I�ɐݒ肳����
*         ���͂̊�_��݂��Ȃ��ꍇ (IRELP=0,1)
*         ���L(1)��������(2)�̂����ꂩ��I��
*         ���͂̊�_��݂��Ȃ��ꍇ (IRELP=0,1)
          IF (IRELP.EQ.0. OR. IRELP.EQ.1) PD(IX,IY,IZ) = X(K)
*         1���]���ȉ��̂�����1�����߂���C���͊��݂���ꍇ (IRELP=2)
*         P'(1,1,1)=0 ---> P(1,1,1)=0 �ɃZ�b�g����ꍇ
          IF (IRELP.EQ.2) PD(IX,IY,IZ) = X(K)-X(1)
   90     CONTINUE
   80   CONTINUE
   70 CONTINUE
*���͕␳�̋��E�����̏���
      CALL PDBND
*���x�̏C��
      DO 100 IX = 1,NX-1
        DO 110 IY = 1,NY
          DO 120 IZ = 1,NZ
          UN(IX,IY,IZ)=UN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX+1,IY,IZ))/DX*DT
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
      DO 130 IX = 1,NX
        DO 140 IY = 1,NY-1
          DO 150 IZ = 1,NZ
          VN(IX,IY,IZ)=VN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX,IY+1,IZ))/DY*DT
  150     CONTINUE
  140   CONTINUE
  130 CONTINUE
      DO 160 IX = 1,NX
        DO 170 IY = 1,NY-1
          DO 180 IZ = 1,NZ
          WN(IX,IY,IZ)=WN(IX,IY,IZ)+(PD(IX,IY,IZ)-PD(IX,IY,IZ+1))/DZ*DT
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*�V���ɓ���ꂽ���x��p���ċ��E��������������
      CALL VELBND
*���͂̏C��
      DO 190 IX = 1,NX
        DO 200 IY = 1,NY
          DO 210 IZ = 1,NZ
            PN(IX,IY,IZ) = PO(IX,IY,IZ) + PD(IX,IY,IZ)
  210     CONTINUE
  200   CONTINUE
  190 CONTINUE
*�V���ɓ���ꂽ���͂�p���ċ��E��������������
      CALL PBND
      RETURN
      END
*********************************************************************
*                     ���x��̌v�Z
*********************************************************************
      SUBROUTINE CALTEM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     ���͕␳�̐��`�V�X�e���p
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     ���͕␳�̐��`�V�X�e���p
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*     ���͕␳�̐��`�V�X�e���p
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0),X(NE0)
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
*                                                                   *
*                ���͕␳�̌W���s��Ƌ��E����                       *
*                                                                   *
* �{�v�Z�v���O�����ł́C���͕␳�̋��E�������C�W���s��쐬���ɍl��  *
* ����D�܂��C����1���Ɨ����̊m�ۂ������ŏ�������D                 *
*                                                                   *
*********************************************************************
      SUBROUTINE PDBNDC (A1,A2,A3,A4,A5,A6,A7,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
      COMMON / D1 / NX,NY,NZ
      COMMON / D2 / DX,DY,DZ,DT
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      DIMENSION A1(NE0),A2(NE0),A3(NE0),A4(NE0),A5(NE0),
     $          A6(NE0),A7(NE0),B(NE0)
*
C     --- SMAC ---
*     �W���s��ɋ��E�����𔽉f������
*     �v�Z�̈�̍����̌W���s��i���E�����𔽉f�j
      IX = 1
      DO 10 IY = 2,NY-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 20 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   20   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   10 CONTINUE
*     �v�Z�̈�̉E���̌W���s��i���E�����𔽉f�j
      IX = NX
      DO 30 IY = 2,NY-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 40 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   40   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   30 CONTINUE
*     �v�Z�̈�̉����̌W���s��i���E�����𔽉f�j
      IY = 1
      DO 50 IX = 2,NX-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 60 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   60   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   50 CONTINUE
*     �v�Z�̈�̏㑤�̌W���s��i���E�����𔽉f�j
      IY = NY
      DO 70 IX = 2,NX-1
        IZ=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 80 IZ = 2,NZ-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   80   CONTINUE
        IZ=NZ
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   70 CONTINUE
*     �v�Z�̈�̌�ʂ̌W���s��i���E�����𔽉f�j
      IZ = 1
      DO 90 IX = 2,NX-1
        IY=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 100 IY = 2,NY-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  100   CONTINUE
        IY=NY
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A7(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
   90 CONTINUE
*     �v�Z�̈�̑O�ʂ̌W���s��i���E�����𔽉f�j
      IZ = NZ
      DO 110 IX = 2,NX-1
        IY=1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
        DO 120 IY = 2,NY-1
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
          A4(I) = - 1.0D0 / DY**2 * DT
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  120   CONTINUE
        IY=NY
          I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
          A3(I) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
          A2(I) = - 1.0D0 / DY**2 * DT
          A1(I) = - 1.0D0 / DX**2 * DT
          A5(I) = - 1.0D0 / DX**2 * DT
          A6(I) = - 1.0D0 / DZ**2 * DT
          B(I) = -DIV(IX,IY,IZ)
  110 CONTINUE
*     �v�Z�̈�̍���_��̌W���s��i���E�����𔽉f�j
      IX=1
      IY=NY
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 130 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  130 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     �v�Z�̈�̍����_��̌W���s��i���E�����𔽉f�j
      IX=1
      IY=1
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 140 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  140 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A5(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     �v�Z�̈�̉E��_��̌W���s��i���E�����𔽉f�j
      IX=NX
      IY=NY
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 150 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  150 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A2(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*     �v�Z�̈�̉E���_��̌W���s��i���E�����𔽉f�j
      IX=NX
      IY=1
      IZ=1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
      DO 160 IZ = 2,NZ-1
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 2.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        A7(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
  160 CONTINUE
      IZ=NZ
        I = IY + (IX-1)*NY + (IZ-1)*(NX*NY)
        A3(I) = DT*( 1.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(I) = - 1.0D0 / DY**2 * DT
        A1(I) = - 1.0D0 / DX**2 * DT
        A6(I) = - 1.0D0 / DZ**2 * DT
        B(I) = -DIV(IX,IY,IZ)
*�P���Ɨ��ȉ��𓾂邽�߂̏��� : ���ږ@�ɂ����Ă͕K�{
*(IX=1,IY=1,IZ=1 ---> M=1����_�Ƃ��C���PN(1,1,1)=PD(1,1,1)=0�Ƃ���)
      IF (IRELP.EQ.1) THEN
        A3(1) = 1.0D0
        A4(1) = 0.0D0
        A5(1) = 0.0D0
        A7(1) = 0.0D0
        B(1)  = 0.0D0
*       M=2 �̓_�̏���(M=1�Ƃ̃����N��f��)
        A3(2) = DT*( 1.0D0/DX**2 + 2.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(2) = - 1.0D0 / DY**2 * DT
        A2(2) = 0.0D0
        A5(2) = - 1.0D0 / DX**2 * DT
        A7(2) = - 1.0D0 / DZ**2 * DT
*       M=1+NY �̓_�̏���(M=1�Ƃ̃����N��f��)
        A3(1+NY) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(1+NY) = - 1.0D0 / DY**2 * DT
        A1(1+NY) = 0.0D0
        A5(1+NY) = - 1.0D0 / DX**2 * DT
        A7(1+NY) = - 1.0D0 / DZ**2 * DT
*       M=1+NX*NY �̓_�̏���(M=1�Ƃ̃����N��f��)
        A3(1+NX*NY) = DT*( 2.0D0/DX**2 + 1.0D0/DY**2 + 1.0D0/DZ**2 )
        A4(1+NX*NY) = - 1.0D0 / DY**2 * DT
        A5(1+NX*NY) = - 1.0D0 / DX**2 * DT
        A6(1+NX*NY) = 0.0D0
        A7(1+NX*NY) = - 1.0D0 / DZ**2 * DT
      END IF
*
      RETURN
      END
*********************************************************************
*                                                                   *
*                  ���͕␳�̋��E�����̏���                         *
*                                                                   *
* �{�v�Z�v���O�����ł́C�W���s��쐬���ɋ��E�������l�����Ă���̂ŁC*
* �����I�ɂ́C�ŏI�I�ɉ��z�Z���̒l��z��Ɋi�[���Ă���ɂ����Ȃ��D  *
* ���̃��[�`���͂Ȃ��Ƃ��悢�D                                      *
* �����ߒ��ŋ��E�������l������ꍇ�͕K�{�ŁC�����̂��тɂ��̏�����  *
* �s���K�v��������D                                                *
*                                                                   *
*********************************************************************
      SUBROUTINE PDBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*�E��
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          PD(NX+1,IY,IZ) = PD(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*����
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          PD(0,IY,IZ) = PD(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*���
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          PD(IX,NY+1,IZ) = PD(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*����
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          PD(IX,0,IZ) = PD(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*���
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          PD(IX,IY,NZ+1) = PD(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*�O��
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          PD(IX,IY,0) = PD(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                     ���͂̋��E�����̏���
*********************************************************************
      SUBROUTINE PBND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*�E��
      DO 10 IY = 0,NY+1
        DO 20 IZ = 0,NZ+1
          PN(NX+1,IY,IZ) = PN(NX,IY,IZ)
   20   CONTINUE
   10 CONTINUE
*����
      DO 30 IY = 0,NY+1
        DO 40 IZ = 0,NZ+1
          PN(0,IY,IZ) = PN(1,IY,IZ)
   40   CONTINUE
   30 CONTINUE
*���
      DO 50 IX = 1,NX
        DO 60 IZ = 1,NZ
          PN(IX,NY+1,IZ) = PN(IX,NY,IZ)
   60   CONTINUE
   50 CONTINUE
*����
      DO 70 IX = 1,NX
        DO 80 IZ = 1,NZ
          PN(IX,0,IZ) = PN(IX,1,IZ)
   80   CONTINUE
   70 CONTINUE
*���
      DO 90 IX = 1,NX
        DO 100 IY = 0,NY+1
          PN(IX,IY,NZ+1) = PN(IX,IY,NZ)
  100   CONTINUE
   90 CONTINUE
*�O��
      DO 110 IX = 1,NX
        DO 120 IY = 0,NY+1
          PN(IX,IY,0) = PN(IX,IY,1)
  120   CONTINUE
  110 CONTINUE
      RETURN
      END
*********************************************************************
*                        �f�[�^�o��
*********************************************************************
      SUBROUTINE PROUT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
C     --- SMAC ---
      COMMON/ARRAY2/DIV(NX0,NY0,NZ0),PD(0:NX0+1,0:NY0+1,0:NZ0+1)
      COMMON / D8 / NE
      COMMON / D9 / NITR
*
       WRITE (11) UN
       WRITE (12) VN
       WRITE (13) WN
       WRITE (14) PN,PD
       WRITE (15) TN
      RETURN
      END
*********************************************************************
*                  Tecplot�p�f�[�^�o��
*********************************************************************
      SUBROUTINE TECPLT (FNAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER ( NX0=20, NY0=20, NZ0=20, NE0=8000, NXY0=400 )
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
     $ PO(0:NX0+1,0:NY0+1,0:NZ0+1),PN(0:NX0+1,0:NY0+1,0:NZ0+1),
     $ TO(0:NX0+1,0:NY0+1,0:NZ0+1),TN(0:NX0+1,0:NY0+1,0:NZ0+1)
      CHARACTER FNAME*20
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
*********************************************************************
*                    �e��̐��`�V�X�e����@                         *
*                                                                   *
*  1. �K�E�X�̏����@                                                *
*  2. point-SOR �@                                                  *
*  3. line-SOR �@                                                   *
*  4. �����c���@                                                    *
*  5. BiCGSTAB�@                                                    *
*                                                                   *
*  ��������C3�����̃|�A�\����������7�_�����ߎ��ɂė��U������       *
*  ���`�V�X�e�����������߂̂��̂ŁC�œK�����Ă���                   *
*  ������̃T�u���[�`������������Ƃ��Ă���                         *
*                                                                   *
* [����]                                                            *
*  1. ������������͓K�X�ύX�̂��ƁD                                *
*  2. ������ NX,NY,NE �ƁC�z��錾���� A1,A2,A3,A4,A5,A6,A7,B,X �́C*
*     ����������̃T�u���[�`�����R�[������Ă��� PRESS �ɂ����đΉ� *
*     ������̂Ɠ������O�Ƃ��Ă��邪�CCOMMON���Œ�`���Ă��Ȃ��̂ŁC*
*     �v�Z�@�̒��ł͈قȂ�ϐ��Ƃ��Ē�`�����D�{�v�Z�v���O������  *
*     �����ẮC�ł��邾�����`�V�X�e����@�̃T�u���[�`���ɔėp����  *
*     �������邽�߁C�����āCCOMMON���͎g�p���Ă��Ȃ��D�܂��C������  *
*     �₷�����邽�߁C�T�u���[�`�����R�[������Ă�����Ɠ������O��*
*     ���ꂼ��̈������`���Ă���D�ȍ~���l�D                      *
*                                                                   *
*********************************************************************
*
*********************************************************************
*                                                                   *
*   3�������v���V�A�����U���ɂ��7�_�����ߎ�                        *
*   �ɂē���ꂽ�K���I��Ώ̍s��A���܂񂾐��`�V�X�e��               *
*                        AX=B                                       *
*   ��Gauss�̏����@��p���ĉ����T�u���[�`���D(���I��)             *
*   A�̓o���h�}�g���b�N�X                                           *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*            ---> A(-NXY:NXY,NE)�Ɋi�[���Ȃ���                      *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NZ : z�����i�q������                                           *
*    NXY: NX * NY                                                   *
*    NE : ���i�q�_�� = NX * NY * NZ                                 *
*                                                                   *
*********************************************************************
      SUBROUTINE GB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY,
     $               A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
*
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
      DIMENSION A(-NXY:NXY,NE)
*
*���ږ@�̂Ƃ���(�W���s�񂪓��قłȂ����)�����Ȃ��ŁC�K�����𓾂�
      IFLG = 0
*�}�g���b�N�XA�̃[���N���A
      DO 10 INE = 1,NE
        DO 20 I = -NXY,NXY
          A(I,INE) = 0.0D0
   20   CONTINUE
   10 CONTINUE
*
*�K�v�ȂƂ����A1����A7�܂ł��i�[����
      DO 30 INE = 1,NE
        A(- NY,INE) = A1(INE)
        A(  -1,INE) = A2(INE)
        A(   0,INE) = A3(INE)
        A(   1,INE) = A4(INE)
        A(  NY,INE) = A5(INE)
        A(-NXY,INE) = A6(INE)
        A( NXY,INE) = A7(INE)
   30 CONTINUE
*
*�O�i����
      DO 40 I = 1,NE-1
        IF ( I.LE.NE-NXY ) THEN
          DO 50 J = 1,NXY
            AA = A(-J,I+J)/A(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 60 K = -J+1,NXY-J
              A(K,I+J) = A(K,I+J)-A(N,I)*AA
              N = N + 1
   60       CONTINUE
   50     CONTINUE
        ELSE
          DO 70 J = 1,NE-I
            AA = A(-J,I+J)/A(0,I)
            B(I+J) = B(I+J) - B(I)*AA
            N = 1
            DO 80 K = -J+1,NE-I-J
              A(K,I+J) = A(K,I+J)-A(N,I)*AA
              N = N + 1
   80       CONTINUE
   70     CONTINUE
        END IF
   40 CONTINUE
*
*��ޑ��
*     �W���s��̓��ِ��𔻒�
      IF ( DABS(A(0,NE)).LE.1.0D-50 ) THEN
        WRITE (6,*) ' Matrix is singular : |A(0,NE)| < 1E-50 '
        IFLG = 1
      END IF
      X(NE) = B(NE) / A(0,NE)
*
      DO 90 I = NE-1,1,-1
        S = 0.0D0
        IF ( I.GT.NE-NXY ) THEN
          DO 100 N = 1,NE-I
            S = S + A(N,I)*X(I+N)
  100     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        ELSE
          DO 110 N = 1,NXY
            S = S + A(N,I)*X(I+N)
  110     CONTINUE
          X(I) = ( B(I)-S ) / A(0,I)
        END IF
   90 CONTINUE
*
*�T�u���[�`���I��
      RETURN
      END
*********************************************************************
*                                                                   *
*  point-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*  3�������v���V�A�����U���ɂ��7�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NZ : z�����i�q������                                           *
*    NXY: NX * NY                                                   *
*    NE : ���i�q�_�� = NX * NY * NZ                                 *
*    NITR : ���e������(in3d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in3d.mac)�ɂĐݒ�                *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
*********************************************************************
      SUBROUTINE PSORB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
*
*     IFLG�̏����l��"��������"
      IFLG=1
*
      DO 10 J = 1,NITR
*
        RNORM1 = 0.0D0
*
        I=1
          XOLD = X(I)
          SUM =                           +A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
*
        DO 20 I=2,NY
          XOLD = X(I)
          SUM =               A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   20   CONTINUE
*
        DO 30 I=NY+1,NXY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $                       +A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   30   CONTINUE
*
        DO 40 I=NXY+1,NE-NXY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A6(I)*X(I-NXY)+A7(I)*X(I+NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   40   CONTINUE
*
        DO 50 I=NE-NXY+1,NE-NY
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)+A5(I)*X(I+NY)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   50   CONTINUE
*
        DO 60 I=NE-NY+1,NE-1
          XOLD = X(I)
          SUM = A1(I)*X(I-NY   )+A2(I)*X(I-1)+A4(I)*X(I+1)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
   60   CONTINUE
*
        I=NE
          XOLD = X(I)
          SUM = A1(I)*X(I-NY)+A2(I)*X(I-1)
     $         +A6(I)*X(I-NXY)
          XNEW = ( B(I)-SUM )/A3(I)
          X(I) = XOLD + OMG * ( XNEW - XOLD )
          RNORM1 = RNORM1 + ( XNEW - XOLD )**2
*
*       ��������; �����Ȃ� IFLG=0 �ɐݒ�
        IF (RNORM1.LE.EPSP) THEN
          IFLG=0
          ITR = J
          GO TO 700
        END IF
   10 CONTINUE
*
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  700 CONTINUE
*
*�T�u���[�`���I��
      RETURN
      END
*********************************************************************
*                                                                   *
*  line-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`�� *
*  3�������v���V�A�����U���ɂ��7�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),          *
*                 AT6(NE),AT7(NE)                                   *
*                                                                   *
*    B -> BX(NE) : ���m�x�N�g��                                     *
*    X -> XN(NE) : ���m�x�N�g�� ---> ��������߂�                   *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NZ : z�����i�q������                                           *
*    NXY: NX * NY                                                   *
*    NE : ���i�q�_�� = NX * NY * NZ                                 *
*    NITR : ���e������(in3d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in3d.mac)�ɂĐݒ�                *
*           1.0�ŏ\���D                                             *
*           ���ӁFPoint-SOR�ƈقȂ�C���܂�傫����������Ɣ��U���� *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
* [�z��̐���]                                                      *
* XN...�e�����ւ̑|�����X(�ԍ��t���͕s��)                          *
*      �͂��߂ɂ��̃T�u���[�`���֓n�����X�ł�����                  *
* X1...�e�����ւ̑|�����X(�ԍ��t���͎������ɈقȂ�)                *
* XO...�e�����ւ̑|���O��X(�ԍ��t���͕s��)                          *
* XOLD...���̃T�u���[�`���ɓ���O��X(�ԍ��t���͕s��)                *
*                                                                   *
*********************************************************************
      SUBROUTINE LSORB (AT1,AT2,AT3,AT4,AT5,AT6,AT7,
     $                  BX,XN,NX,NY,NZ,NE,NXY,
     $                  X1,XO,XOLD,A,B,C,D,U,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION AT1(NE),AT2(NE),AT3(NE),AT4(NE),AT5(NE),
     $          AT6(NE),AT7(NE)
      DIMENSION XN(NE),X1(NE),BX(NE),XO(NE),XOLD(NE)
      DIMENSION A(NE),B(NE),C(NE),D(NE),U(NE),Y(NE)
*
*     IFLG�̏����l��"��������"
      IFLG=1
*
      DO 10 K=1,NITR
*     x �������ւ̑|�� : �g�[�}�X�@�ɂ��
      INX = 1
      DO 100 IZ = 1,NZ
        DO 110 IY = 1,NY
          DO 120 IX = 1,NX
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X
          A(INX) = AT1(INZ)
          B(INX) = AT3(INZ)
          C(INX) = AT5(INZ)
          D(INX) = BX(INZ)
          IF (INZ-1.GE.1) THEN
            D(INX)=D(INX)-AT2(INZ)*XN(INZ-1)
          END IF
          IF (INZ+1.LE.NE) THEN
            D(INX)=D(INX)-AT4(INZ)*XN(INZ+1)
          END IF
          IF (INZ-NXY.GE.1) THEN
            D(INX)=D(INX)-AT6(INZ)*XN(INZ-NXY)
          END IF
          IF (INZ+NXY.LE.NE) THEN
            D(INX)=D(INX)-AT7(INZ)*XN(INZ+NXY)
          END IF
*         �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ�
          XO(INZ) = XN(INZ)
*         x�����ւ̑|���̑O�̒l��XOLD�ɕۑ�
          XOLD(INZ) = XN(INZ)
          INX = INX + 1
  120     CONTINUE
  110   CONTINUE
  100 CONTINUE
*     Ly=b ������
      U(1) = C(1) / B(1)
      DO 130 J = 2,NE-1
        U(J) = C(J) / ( B(J)-A(J)*U(J-1) )
  130 CONTINUE
      Y(1) = D(1) / B(1)
      DO 140 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  140 CONTINUE
*     Ux=y ������
      X1(NE) = Y(NE)
      DO 150 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  150 CONTINUE
      INX = 1
      DO 160 IZ = 1,NZ
        DO 170 IY = 1,NY
          DO 180 IX = 1,NX
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INX)
            INX = INX + 1
  180     CONTINUE
  170   CONTINUE
  160 CONTINUE
*     y �������ւ̑|�� : �g�[�}�X�@�ɂ��
      INY = 1
      DO 200 IZ = 1,NZ
        DO 210 IX = 1,NX
          DO 220 IY = 1,NY
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X
          A(INY) = AT2(INZ)
          B(INY) = AT3(INZ)
          C(INY) = AT4(INZ)
          D(INY) = BX(INZ)
          IF (INZ-NY.GE.1) THEN
            D(INY)=D(INY)-AT1(INZ)*XN(INZ-NY)
          END IF
          IF (INZ+NY.LE.NE) THEN
            D(INY)=D(INY)-AT5(INZ)*XN(INZ+NY)
          END IF
          IF (INZ-NXY.GE.1) THEN
            D(INY)=D(INY)-AT6(INZ)*XN(INZ-NXY)
          END IF
          IF (INZ+NXY.LE.NE) THEN
            D(INY)=D(INY)-AT7(INZ)*XN(INZ+NXY)
          END IF
*         �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ�
          XO(INZ) = XN(INZ)
          INY = INY + 1
  220     CONTINUE
  210   CONTINUE
  200 CONTINUE
*     Ly=b ������
      U(1) = C(1) / B(1)
      DO 230 J = 2,NE-1
        U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  230 CONTINUE
      Y(1) = D(1) / B(1)
      DO 240 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  240 CONTINUE
*     Ux=y ������
      X1(NE) = Y(NE)
      DO 250 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  250 CONTINUE
      INY = 1
      DO 260 IZ = 1,NZ
        DO 270 IX = 1,NX
          DO 280 IY = 1,NY
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INY)
            INY = INY + 1
  280     CONTINUE
  270   CONTINUE
  260 CONTINUE
*
*     z �������ւ̑|�� : �g�[�}�X�@�ɂ��
      INZZ = 1
      DO 300 IY = 1,NY
        DO 310 IX = 1,NX
          DO 320 IZ = 1,NZ
          INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*         �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X
          A(INZZ) = AT6(INZ)
          B(INZZ) = AT3(INZ)
          C(INZZ) = AT7(INZ)
          D(INZZ) = BX(INZ)
          IF (INZ-NY.GE.1) THEN
            D(INZZ)=D(INZZ)-AT1(INZ)*XN(INZ-NY)
          END IF
          IF (INZ+NY.LE.NE) THEN
            D(INZZ)=D(INZZ)-AT5(INZ)*XN(INZ+NY)
          END IF
          IF (INZ-1.GE.1) THEN
            D(INZZ)=D(INZZ)-AT2(INZ)*XN(INZ-1)
          END IF
          IF (INZ+1.LE.NE) THEN
            D(INZZ)=D(INZZ)-AT4(INZ)*XN(INZ+1)
          END IF
*         �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ�
          XO(INZ) = XN(INZ)
          INZZ = INZZ + 1
  320     CONTINUE
  310   CONTINUE
  300 CONTINUE
*     Ly=b ������
      U(1) = C(1) / B(1)
      DO 330 J = 2,NE-1
        U(J) = C(J)/( B(J)-A(J)*U(J-1) )
  330 CONTINUE
      Y(1) = D(1) / B(1)
      DO 340 J = 2,NE
        Y(J) = ( D(J)-A(J)*Y(J-1) ) / ( B(J)-A(J)*U(J-1) )
  340 CONTINUE
*     Ux=y ������
      X1(NE) = Y(NE)
      DO 350 J = NE-1,1,-1
        X1(J) = Y(J) - U(J)*X1(J+1)
  350 CONTINUE
      INZZ = 1
      DO 360 IY = 1,NY
        DO 370 IX = 1,NX
          DO 380 IZ = 1,NZ
            INZ = IY + (IX-1)*NY + (IZ-1)*NXY
*           ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa
            XN(INZ)=(1.0D0-OMG)*XO(INZ)+OMG*X1(INZZ)
            INZZ = INZZ + 1
  380     CONTINUE
  370   CONTINUE
  360 CONTINUE
*
      RNORM= 0.0D0
      DO 400 I = 1,NE
        RNORM = RNORM + (XN(I)-XOLD(I))**2
  400 CONTINUE
*
*     ��������
      IF (RNORM.LE.EPSP) THEN
        IFLG=0
        ITR=K
        GO TO 900
      END IF
*
   10 CONTINUE
*
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  900 CONTINUE
*
      RETURN
      END
*********************************************************************
*                                                                   *
*    �����c��(Conjugate Residual)�@�ɂ���Ώ̍s�� A ���܂�        *
*  ���`�V�X�e����@�T�u���[�`��                                     *
*                        AX=B                                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*                                                                   *
*  B : ���m�x�N�g��                                                 *
*  X : ���m�x�N�g�� ---> ��������߂� ---> �����ł͕֋X��z�� XP(NE)*
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NZ : z�����i�q������                                           *
*    NXY: NX * NY                                                   *
*    NE : ���i�q�_�� = NX * NY * NZ                                 *
*    NITR : ���e������(in3d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in3d.mac)�ɂĐݒ�                *
*                                                                   *
*�m�z��̐����n                                                     *
*    R(NE) : r_{k} = B - A x_{k}                                    *
*    P(NE) :  p_{k+1} = r_{k+1} + ��_{k} p_{k}, p_{0} = r_{0}       *
*    AP(NE) : A * P                                                 *
*    AR(NE) : A * R                                                 *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
*********************************************************************
      SUBROUTINE CRB (A1,A2,A3,A4,A5,A6,A7,B,XP,NX,NY,NZ,NE,NXY,
     $                R,P,AP,AR,X,XOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),XP(NE)
*��Ɨp�z��
      DIMENSION R(NE),P(NE),AP(NE),AR(NE),X(NE),XOLD(NE)
*
* R �� AX ����
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,X,R,NX,NY,NZ,NE,NXY)
*
      DO 40 I = 1, NE
*       r_{0}��p_{0}(�����l)�̐ݒ�
        R(I) = B(I) - R(I)
        P(I) = R(I)
*       �O�̎�����X��XOLD�ɑ��
        XOLD(I) = XP(I)
   40 CONTINUE
*
* AP�� A p_{0} ����
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NX,NY,NZ,NE,NXY)
*
* �����v�Z
      DO 50 K = 1,NITR
*       ( r_{k}, A p_{k} )�̌v�Z => RAP
        CALL PROVV (R,AP,RAP,NE)
*       ( A p_{k}, A p_{k} )�̌v�Z => APAP
        CALL PROVV (AP,AP,APAP,NE)
*       ��_{k} = ( r_{k}, A p_{k} ) / ( A p_{k}, A p_{k} )
**********************************************************
*.......�T�������̂��߂̌v�Z�� 0���Z �Ȃ�v�Z�I��(--- SMAC ---)
        IF (DABS(APAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : ALPHA_{K} in Conjugate Residual '
          IFLG = 2
          RETURN
        ELSE
          ALP = RAP / APAP
        END IF
**********************************************************
*
        RNORM = 0.0D0
*
        DO 70 I = 1,NE
*         x_{k+1}=x_{k}+��_{k}p_{k}
          X(I) = X(I) + ALP*P(I)
*         r_{k+1}=r_{k}-��_{k}Ap_{k}
          R(I) = R(I) - ALP*AP(I)
*         �O�̔����Ƃ̍��̃m�����̌v�Z
          RNORM = RNORM + (X(I)-XOLD(I))**2
*         ����ꂽX��XOLD�ɑ��
          XOLD(I) = X(I)
   70   CONTINUE
*
* RNORM �� EPSP �ȉ��Ȃ�����Ƃ݂Ȃ��� 700 ��
        IF (RNORM.LE.EPSP) THEN
          IFLG=0
          ITR=K
          GO TO 700
        END IF
*
* ���������̏ꍇ
*       A r_{k+1} �̌v�Z => AR(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,R,AR,NX,NY,NZ,NE,NXY)
*       ( A r_{k+1}, A p_{k} )�̌v�Z => ARAP
        CALL PROVV( AR,AP,ARAP,NE)
*       ��_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} )
**********************************************************
*.......�T�������̂��߂̌v�Z�� 0���Z �Ȃ�v�Z�I��(--- SMAC ---)
        IF (DABS(APAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : BETA_{K} in Conjugate Residual'
          IFLG = 2
          RETURN
        ELSE
          BETA = -ARAP / APAP
        END IF
**********************************************************
*
        DO 90 I = 1, NE
*         p_{k+1} = r_{k+1} + ��_{k} p_{k}
          P(I) = R(I) + BETA*P(I)
*         A p_{k+1} = A r_{k+1} + ��_{k}A p_{k}
          AP(I) = AR(I) + BETA*AP(I)
   90   CONTINUE
*
   50 CONTINUE
*     NITR �܂Ōv�Z���Ă���������
      IFLG=1
*
* �����Ɣ��肳�ꂽ�Ƃ��̕���_
  700 CONTINUE
*
      DO 100 I = 1,NE
        XP(I) = X(I)
  100 CONTINUE
*
*�T�u���[�`���I��
      RETURN
      END
*
*********************************************************************
*                                                                   *
*    �x�N�g�� A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��                 *
*                        AB=C                                       *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NE : ���i�q�_��(�x�N�g�� A,B �̃T�C�Y)                         *
*    C  : A �� B �̐�(�X�J���[)                                     *
*                                                                   *
*********************************************************************
      SUBROUTINE PROVV(A,B,C,NE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(NE),B(NE)
*
      C = 0.0D0
      DO 10 I=1,NE
        C = C + A(I)*B(I)
   10 CONTINUE
      RETURN
      END
*
*********************************************************************
*                                                                   *
*    �}�g���b�N�X A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��             *
*                        AB=C                                       *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NE : ���i�q�_��(�����}�g���b�N�X A,B,C �̃T�C�Y)               *
*    C  : A �� B �̐�(�x�N�g��)                                     *
*                                                                   *
*********************************************************************
      SUBROUTINE PROMV (A1,A2,A3,A4,A5,A6,A7,B,C,NX,NY,NZ,NE,NXY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),C(NE)
*
      I=1
        C(I) =                              A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
*
      DO 10 I=2,NY
        C(I) = A2(I)*B(I-1)                +A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
   10 CONTINUE
*
      DO 20 I=NY+1,NXY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $                       +A7(I)*B(I+NXY)
   20 CONTINUE
*
      DO 30 I=NXY+1,NE-NXY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $        +A6(I)*B(I-NXY)+A7(I)*B(I+NXY)
   30 CONTINUE
*
      DO 40 I=NE-NXY+1,NE-NY
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)  +A5(I)*B(I+NY)
     $        +A6(I)*B(I-NXY)
   40 CONTINUE
*
      DO 50 I=NE-NY+1,NE-1
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A4(I)*B(I+1)
     $        +A6(I)*B(I-NXY)
   50 CONTINUE
*
      I=NE
        C(I) = A1(I)*B(I-NY) +A2(I)*B(I-1 )+A3(I)*B(I)
     $        +A6(I)*B(I-NXY)
*
*�T�u���[�`���I��
      RETURN
      END
*
*********************************************************************
*                                                                   *
*  Bi-CGSTAB �@�ɂ���Ώ̍s�� A ���܂�                            *
*  ���`�V�X�e����@�T�u���[�`��                                     *
*                        AX=B                                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE)  *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NZ : z�����i�q������                                           *
*    NXY: NX * NY                                                   *
*    NE : ���i�q�_�� = NX * NY * NZ                                 *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in3d.mac)�ɂĐݒ�                *
*                                                                   *
*�m�z��̐����n                                                     *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
*********************************************************************
      SUBROUTINE BICGB (A1,A2,A3,A4,A5,A6,A7,B,X,NX,NY,NZ,NE,NXY,
     $                  R,AT,AP,P,S,T,XOLD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON / D4 / RE,PR,GR,TIME,OMG,EPSP
      COMMON / D5 / ICYCLE,ITR,IFLG,IRELP,METHOD
C     --- SMAC ---
      COMMON / D9 / NITR
*
      DIMENSION A1(NE),A2(NE),A3(NE),A4(NE),A5(NE),A6(NE),A7(NE),
     $          B(NE),X(NE)
* ��Ɨp�z��
      DIMENSION R(NE),AT(NE),AP(NE),P(NE),S(NE),T(NE),XOLD(NE)
*
      DO 5 I = 1,NE
        XOLD(I) = X(I)
    5 CONTINUE
*
* R �� AX ����
      CALL PROMV (A1,A2,A3,A4,A5,A6,A7,X,R,NX,NY,NZ,NE,NXY)
* r_{0}��p_{0}(�����l)�C������ s=r_{0} �̐ݒ�
      DO 10 I=1,NE
        R(I) = B(I) - R(I)
        P(I) = R(I)
        S(I) = R(I)
   10 CONTINUE
*
* �J��Ԃ��v�Z
      DO 20 J =1,NITR
*       ( s, r_{k} ) �̌v�Z => SR1
        CALL PROVV (S,R,SR1,NE)
*       A p_{k} �̌v�Z => AP(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,P,AP,NX,NY,NZ,NE,NXY)
*       ( s, A p_{k} ) �̌v�Z => SAP
        CALL PROVV (S,AP,SAP,NE)
*       ��_{k} = ( s, r_{k} ) / ( s, A p_{k} )
**********************************************************
*.......�T�������̂��߂̌v�Z�� 0���Z �Ȃ�v�Z�I��(--- SMAC ---)
        IF (DABS(SAP).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : ALPHA_{K} in Bi-CGSTAB'
          IFLG = 2
          RETURN
        ELSE
          ALPHA = SR1/SAP
        END IF
**********************************************************
        DO 50 I=1,NE
*         t_{k} = r_{k} - ��_{k} A p_{k}
          T(I) = R(I) - ALPHA*AP(I)
   50   CONTINUE
*       A t_{k} �̌v�Z => AT(NE)
        CALL PROMV (A1,A2,A3,A4,A5,A6,A7,T,AT,NX,NY,NZ,NE,NXY)
*       ( A t_{k}, t_{k} ) �̌v�Z => ATT
        CALL PROVV (AT,T,ATT,NE)
*       ( A t_{k}, A t_{k} ) �̌v�Z => ATAT
        CALL PROVV (AT,AT,ATAT,NE)
*       ��_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} )
**********************************************************
*.......�T�������̂��߂̌v�Z�� 0���Z �Ȃ�v�Z�I��(--- SMAC ---)
        IF (DABS(ATAT).LT.1.0D-50) THEN
          WRITE (6,*) ' 0 division : XI_{K} in Bi-CGSTAB'
          IFLG = 2
          RETURN
        ELSE
          XI = ATT/ATAT
        END IF
**********************************************************
        RNORM = 0.0D0
        DO 60 I=1,NE
*         x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{k} t_{k}
          X(I) = X(I) + ALPHA*P(I) + XI*T(I)
*         r_{k+1} = t_{k} - ��_{k} A t_{k}
          R(I) = T(I) - XI*AT(I)
*         �O�̔����Ƃ̍��̃m�����̌v�Z
          RNORM = RNORM + (X(I)-XOLD(I))**2
*         ����ꂽX��XOLD�ɑ��
          XOLD(I) = X(I)
   60   CONTINUE
* RNORM �� EPSP �ȉ��Ȃ�����Ƃ݂Ȃ��� 900 ��
        IF (RNORM.LE.EPSP) THEN
          IFLG=0
          ITR=J
          GO TO 900
        END IF
*       ���������̏ꍇ
*       ( s, r_{k+1} ) �̌v�Z => SR2
        CALL PROVV (S,R,SR2,NE)
*       ��_{k} = ( ��_{k} / ��_{k} ) * ( s, r_{k+1} )/( s, r_{k} )
        BETA = (ALPHA / XI) * (SR2 / SR1)
        DO 70 I=1,NE
*         p_{k+1} = r_{k+1} + ��_{k}( p_{k}-��_{k} A p_{k} )
          P(I) = R(I) + BETA * ( P(I) - XI*AP(I) )
   70   CONTINUE
   20 CONTINUE
*     NITR �܂Ōv�Z���Ă���������
      IFLG=1
*
  900 CONTINUE
      RETURN
      END
*
