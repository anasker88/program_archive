/********************************************************************
* �t�@�C����  �Fsmac3d.c                                            *
* �^�C�g��    �FSMAC�@�ɂ��3�����M������̓v���O����               *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� �o�C�I�E���p���w��              *
* �����      �F2011.12.25                                          *
* ����        �FC                                                   *
*********************************************************************
* �v���O���������s����ƁC"in3d.mac"(�K���������j��ǂ݂ɍs���D     *
* �ڍׂ�FORTRAN�v���O�����̃R�����g���Q�ƁD                         *
*                                                                   *
*  ���ӁF�b����ł͔z��̓Y�����O����n�܂�̂ŁCFORTRAN ��X(N)��   *
*  ���̂܂�X[N]�Ɛ錾����ƁCX[0]����X[N-1]���Ӗ�����D�����ł�     *
*  �ł������FORTRAN�v���O�����ɋ߂��Ȃ�悤�ɁCX[N+1]�Ƃ��ēY��    *
*  ��1����N�܂ŕω������Ă���D                                     *
********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NX 20 /* NX:x�����i�q�� */
#define NY 20 /* NY:y�����i�q�� */
#define NZ 20 /* NZ:z�����i�q�� */
#define NE 8000 /* NE:���͕␳�̖��m���̐� --- SMAC --- */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NZ1 21 /* NZ1=NZ+1 */
#define NZ2 22 /* NZ2=NZ+2 */
#define NE1 8001 /* NE1=NE+1 */
#define NXY 400 /* NXY=NX*NY --- SMAC --- */
#define NXYNXY 801 /* NXYNXY=2*NXY+1 �K�E�X�̏����@�Ŏg�p --- SMAC --- */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
void pdbndc (double [NE1], double [NE1], double [NE1], double [NE1], double [NE1],
             double [NE1], double [NE1], double [NE1]);
void gb (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
         double [NE1], double [NE1], double [NE1], double [NE1],
           int nx, int ny, int nz, int ne, int nxy);
void crb   (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
            double [NE1], double [NE1], double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void promv (double [NE1], double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void bicgb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1],
            double [NE1], double [NE1], 
            int nx, int ny, int nz, int ne, int nxy);
double provv (double [NE1], double [NE1], int ne);
void psorb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void lsorb (double [NE1] ,double [NE1], double [NE1], double [NE1], double [NE1], 
            double [NE1], double [NE1],
            double [NE1], double [NE1],
            int nx, int ny, int nz, int ne, int nxy);
void pdbnd();
void pbnd();
double DX, DY, DZ, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
int NITR;
char fname[12][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_15, *out_21;
FILE *in_16, *in_17, *in_18, *in_19, *in_20;
double UO[NX1][NY2][NZ2],UN[NX1][NY2][NZ2],VO[NX2][NY1][NZ2],VN[NX2][NY1][NZ2];
double WO[NX2][NY2][NZ1],WN[NX2][NY2][NZ1],PO[NX2][NY2][NZ2],PN[NX2][NY2][NZ2];
double TO[NX2][NY2][NZ2],TN[NX2][NY2][NZ2];
double DIV[NX1][NY1][NZ1], PD[NX2][NY2][NZ2];/* --- SMAC ---*/
void main()
{
     int i,j,NCYCLE;
     double DLX,DLY,DLZ;
     char buff[80];
     fname10="in3d.mac";
     strncpy(fname[0],fname10,8);
     in_10=fopen(fname[0],"rt"); /*�p�����[�^�t�@�C���̃I�[�v��*/
     for (i=1;i<=11;i++){/*�o�̓t�@�C�����̓ǂݍ���*/
       fgets(buff, sizeof buff,in_10);
       for (j=0;buff[j]!=0;j++);
       strncpy(fname[i], buff, j-1);
     }
     out_11=fopen(fname[1],"wb");/*U�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_12=fopen(fname[2],"wb");/*V�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_13=fopen(fname[3],"wb");/*W�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_14=fopen(fname[4],"wb");/*P�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_15=fopen(fname[5],"wb");/*T�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(12�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*                   ��(13�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(14�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(15�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(16�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(17�s��)*/
     fscanf(in_10," %d %d %d %d ",&ITYPE,&ICYCLE,&NITR,&NCYCLE);
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(19�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(20�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg ",&EPSP,&OMG);
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(22�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(23�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg %lg %lg ",&DT,&RE,&PR,&GR);
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(25�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in3d.mac���̃R�����g�s(26�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg %lg %d %d",&DLX,&DLY,&DLZ,&IRELP,&METHOD);
     printf (" ITYPE = %d ICYCLE= %d NITR= %d NCYCLE= %d \n",
             ITYPE,ICYCLE,NITR,NCYCLE);
     printf(" EPSP= %12.3lE  OMG = %12.3lE \n",EPSP,OMG);
     printf(" DT = %12.3lE  RE = %12.3lE  PR = %12.3lE  GR = %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX=%12.3lE DLY=%12.3lE DLZ=%12.3lE IRELP=%d METHOD=%d \n",DLX,DLY,DLZ,IRELP,METHOD);
     if ( ICYCLE != 0 ){/*�p���̌v�Z�̏ꍇ*/
       in_16=fopen(fname[6],"rb");/*U�f�[�^�t�@�C���̃I�[�v��*/
       in_17=fopen(fname[7],"rb");/*V�f�[�^�t�@�C���̃I�[�v��*/
       in_18=fopen(fname[8],"rb");/*W�f�[�^�t�@�C���̃I�[�v��*/
       in_19=fopen(fname[9],"rb");/*P�f�[�^�t�@�C���̃I�[�v��*/
       in_20=fopen(fname[10],"rb");/*T�f�[�^�t�@�C���̃I�[�v��*/
     }
     DX=DLX/(double)NX; DY=DLY/(double)NY; DZ=DLZ/(double)NZ;/*x,y,z�����̊i�q��DX,DY,DZ*/
     VIS=PR; ALP=1.0; BUO=( GR * ( PR*PR ) );
      /*VIS:�^�����������̊g�U���̌W��(�����ł�Pr)*/
      /*ALP:�G�l���M�[���������̊g�U���̌W��(�����ł�1)*/
      /*BUO:���͍��̌W��(�����ł� Gr * Pr**2)*/
     if ( ITYPE==1 ){/* ������Ȃ畂�͍��̌W�����[���ɐݒ� */
       BUO = 0.0;
     }
     cinit();/*�����l�̐ݒ�*/
     label_700:{};/*���Ԑi�s�̂��߂̖߂�_*/
     adv();/*���Ԑi�s*/
     IFLG=0;
     /*  --- SMAC ---
      �{�v�Z�ɂ����ẮCICYCLE=0�̂Ƃ��̏��������͑��x��ƈ��͏���[����
      ���Ă���̂ŁC���͕␳�̐��`�V�X�e���������̂�����ƂȂ邽�߁C
      ICYCLE=1�̂Ƃ������͉��x�̌v�Z�փW�����v����悤�ɂ���D
      �v�Z��������ɉ����ēK�X�폜���邢�͕ύX */
      if ( ICYCLE==1 ){
        goto label_720;
      }
     calvel();/*���x��̌v�Z*/
     press();/*���̈��͏�̌v�Z*/
     label_720:{};/* ���x��ƈ��͏ꂪ�[���ƂȂ�Ƃ��̃X�L�b�v�� --- SMAC --- */
     if ( IFLG==0 ){/*���͏�̌v�Z�����������Ƃ�*/
       if ( ITYPE==2 ){/*�񓙉���v�Z�̏ꍇ*/
         caltem();/*���x����v�Z*/
       }
     }
     if ( IFLG==1 ){/*���͏�̌v�Z���������Ă��Ȃ��Ƃ�*/
       /*�f�[�^���o�͂��ċ����I��*/
       printf(" calculation has diverged \n");
       prout();
       goto label_900;
     }
     if ( ICYCLE < NCYCLE ){/*���Ԑi�s�J�E���^(ICYCLE)��NCYCLE��菬������*/
       goto label_700;
     }
     else{/*���Ԑi�s�J�E���^��NCYCLE�ɂȂ�����->�v�Z�I��*/
       prout();
     }
     label_900:{};
     tecplt(fname[11]);/*Tecplot�p�f�[�^�̏o��*/
     fclose (in_10);  fclose (out_11); fclose (out_12);
     fclose (out_13); fclose (out_14); fclose(out_15);
}
/*�����ݒ�*/
void cinit()
{
     int ix,iy,iz;
     if ( ICYCLE == 0 ){/*�V�K�v�Z�̏ꍇ*/
       for (ix=0;ix<=NX;ix++){/*U�̏����l�ݒ�*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             UN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*V�̏����l�ݒ�*/
         for (iy=0;iy<=NY;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             VN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*W�̏����l�ݒ�*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ;iz++){
             WN[ix][iy][iz] = 0.0;
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*P�̏����l�ݒ�*/
         for (iy=0;iy<=NY+1;iy++){
           for (iz=0;iz<=NZ+1;iz++){
             PD[ix][iy][iz] = 0.0; PN[ix][iy][iz] = 0.0;/* --- SMAC ---*/
           }
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*T�̏����l�ݒ�(�̈���͍���(+0.5)�ƒቷ(-0.5)�̒��ԉ��x)  */
         for (iy=0;iy<=NY+1;iy++){/* �i���Ӂj���͍��̌v�Z�ŉ��x�̔z����g�p���Ă���̂ŁC  */
           for (iz=0;iz<=NZ+1;iz++){/*       ������ł�T=0�Ƃ��ď������������͐ݒ肷��K�v */
             TN[ix][iy][iz] = 0.0;/*         ������D�[���ȊO�̒l������ƕ��͍����v�Z����*/
           }/*                               ��\��������̂Œ��ӁD*/
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*T�̋��E�F�E���ǁi��p�jT=-0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[NX+1][iy][iz] = 2.0 * ( -0.5 ) - TN[NX][iy][iz];
         }
       }
       for (iy=0;iy<=NY+1;iy++){/*T�̋��E�F�����ǁi���M�jT=+0.5*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[0][iy][iz] = 2.0 * ( +0.5 ) - TN[1][iy][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F��ʁi�f�M�j*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[ix][NY+1][iz] = TN[ix][NY][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F���ʁi�f�M�j*/
         for (iz=0;iz<=NZ+1;iz++){
           TN[ix][0][iz] = TN[ix][1][iz];
         }
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F��ʁi�f�M�j*/
         for (iy=1;iy<=NY;iy++){
           TN[ix][iy][NZ+1] = TN[ix][iy][NZ];
         }
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F�O�ʁi�f�M�j*/
         for (iy=1;iy<=NY;iy++){
           TN[ix][iy][0] = TN[ix][iy][1];
         }
       }
     }
     else{/*�p���v�Z�i���łɂ���v�Z���ʂ���X�^�[�g�j�̏ꍇ*/
       /*U�f�[�^�t�@�C������̓ǂݍ���[Unit No.=16]*/
       fread(UN, sizeof(double), NX1*NY2*NZ2, in_16);
       /*V�f�[�^�t�@�C������̓ǂݍ���[Unit No.=17]*/
       fread(VN, sizeof(double), NX2*NY1*NZ2, in_17);
       /*W�f�[�^�t�@�C������̓ǂݍ���[Unit No.=18]*/
       fread(WN, sizeof(double), NX2*NY2*NZ1, in_18);
       /*P�f�[�^�t�@�C������̓ǂݍ���[Unit No.=19]*/
       fread(PN, sizeof(double), NX2*NY2*NZ2, in_19);
       fread(PD, sizeof(double), NX2*NY2*NZ2, in_19);
       /*T�f�[�^�t�@�C������̓ǂݍ���[Unit No.=20]*/
       fread(TN, sizeof(double), NX2*NY2*NZ2, in_20);
       fclose (in_16); fclose (in_17); fclose (in_18); fclose (in_19); fclose (in_20);
     }
}
/*���Ԑi�s*/
void adv()
{
     int ix,iy,iz;
     TIME = DT*(double)ICYCLE; ICYCLE = ICYCLE + 1;
     if ( (ICYCLE%100) ==0 ){/*ICYCLE��100�񖈂ɕ\��*/
       printf ("  CYC = %8d \n",ICYCLE);
     }
     for (ix=0;ix<=NX;ix++){/* UN -> UO : �K�v�Ȃ����ւ���O��UN��UO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/* UN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         for (iz=0;iz<=NZ+1;iz++){/*UO:�V�������ԃX�e�b�v�ł̏����l�CUN��ۑ��D*/
           UO[ix][iy][iz] = UN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*    VN -> VO : �K�v�Ȃ����ւ���O��VN��VO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY;iy++){/*    VN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         for (iz=0;iz<=NZ+1;iz++){/*VO:�V�������ԃX�e�b�v�ł̏����l�CVN��ۑ��D*/
           VO[ix][iy][iz] = VN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*WN -> WO : �K�v�Ȃ����ւ���O��WN��WO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*WN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         for (iz=0;iz<=NZ;iz++){/*WO:�V�������ԃX�e�b�v�ł̏����l�CWN��ۑ��D*/
           WO[ix][iy][iz] = WN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : �K�v�Ȃ����ւ���O��TN��TO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*  TN:�O�̎��ԃX�e�b�v�ł̌v�Z�l*/
         for (iz=0;iz<=NZ+1;iz++){/*TO:�V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D*/
           TO[ix][iy][iz] = TN[ix][iy][iz];
           PO[ix][iy][iz] = PN[ix][iy][iz];/* --- SMAC --- PO:�V�������ԃX�e�b�v�ł̏����l�DPN��ۑ��D*/
         }
       }
     }
}
/*���x��̌v�Z*/
void calvel()
{
     int ix,iy,iz;
     double vv,ww,cnvux,cnvuy,cnvuz,tu,buou,difu;
     double uu,cnvvx,cnvvy,cnvvz,tv,buov,difv;
     double cnvwx,cnvwy,cnvwz,tw,buow,difw;
     /*U(ix,iy)�̌v�Z*/
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* vv��U(ix,iy,iz)�ɂ�����V�̕�Ԓl */
         vv=(VO[ix][iy][iz]+VO[ix+1][iy][iz]+VO[ix][iy-1][iz]+VO[ix+1][iy-1][iz])/4.0;
         /* ww��U(ix,iy,iz)�ɂ�����W�̕�Ԓl */
         ww=(WO[ix][iy][iz]+WO[ix+1][iy][iz]+WO[ix][iy][iz-1]+WO[ix+1][iy][iz-1])/4.0;
         /* �Η���(cnvux,cnvuy)���P�����x���㍷���ɂČv�Z */
         if ( UO[ix][iy][iz] >= 0.0 ){
           cnvux = UO[ix][iy][iz]*( UO[ix][iy][iz]-UO[ix-1][iy][iz] ) / DX;
         }
         else if ( UO[ix][iy][iz] < 0.0 ){
           cnvux = UO[ix][iy][iz]*( UO[ix+1][iy][iz] - UO[ix][iy][iz] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvuy = vv*( UO[ix][iy][iz] - UO[ix][iy-1][iz] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvuy = vv*( UO[ix][iy+1][iz] - UO[ix][iy][iz] ) / DY;
         }
         if ( ww >= 0.0 ){
           cnvuz = ww*( UO[ix][iy][iz] - UO[ix][iy][iz-1] ) / DZ;
         }
         else if ( ww < 0.0 ){
           cnvuz = ww*( UO[ix][iy][iz+1] - UO[ix][iy][iz] ) / DZ;
         }
         /* x�����̕��͍�(buou)�̓[�� */
         tu = 0.0;
         buou = BUO*tu;
         /* �g�U��(difu)�̌v�Z */
         difu = VIS*( (UO[ix-1][iy][iz]-2.0*UO[ix][iy][iz]+UO[ix+1][iy][iz])/(DX*DX)
                     +(UO[ix][iy-1][iz]-2.0*UO[ix][iy][iz]+UO[ix][iy+1][iz])/(DY*DY)
                     +(UO[ix][iy][iz-1]-2.0*UO[ix][iy][iz]+UO[ix][iy][iz+1])/(DZ*DZ) );
         /*���̑��x(U)�̌v�Z*/
         UN[ix][iy][iz] = UO[ix][iy][iz]
              + DT*( -cnvux-cnvuy-cnvuz+difu+buou+( PO[ix][iy][iz]-PO[ix+1][iy][iz] )/DX );
         }
       }
     }
     /*V(ix,iy,iz)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* uu��V(ix,iy,iz)�ɂ�����U�̕�Ԓl */
         uu=(UO[ix-1][iy][iz]+UO[ix][iy][iz]+UO[ix-1][iy+1][iz]+UO[ix][iy+1][iz])/4.0;
         /* ww��V(ix,iy,iz)�ɂ�����W�̕�Ԓl */
         ww=(WO[ix][iy][iz-1]+WO[ix][iy][iz]+WO[ix][iy+1][iz-1]+WO[ix][iy+1][iz])/4.0;
         /* �Η���(cnvvx,cnvvy,cnvvz)���P�����x���㍷���ɂČv�Z */
         if ( uu >= 0.0 ){
           cnvvx = uu*( VO[ix][iy][iz] - VO[ix-1][iy][iz] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvvx = uu*( VO[ix+1][iy][iz] - VO[ix][iy][iz] ) / DX;
         }
         if ( VO[ix][iy][iz] >= 0.0 ){
           cnvvy = VO[ix][iy][iz]*( VO[ix][iy][iz]-VO[ix][iy-1][iz] ) / DY;
         }
         else if ( VO[ix][iy][iz] < 0.0 ){
           cnvvy = VO[ix][iy][iz]*( VO[ix][iy+1][iz]-VO[ix][iy][iz] ) / DY;
         }
         if ( ww >= 0.0 ){
           cnvvz = WO[ix][iy][iz]*( VO[ix][iy][iz]-VO[ix][iy][iz-1] ) / DZ;
         }
         else if ( ww < 0.0 ){
           cnvvz = WO[ix][iy][iz]*( VO[ix][iy][iz+1]-VO[ix][iy][iz] ) / DZ;
         }
         /*���͍�(buov)�̌v�Z*/
         tv = ( TO[ix][iy][iz] + TO[ix][iy+1][iz] )/2.0;
         buov = BUO*tv;
         /*�g�U��(difv)�̌v�Z*/
         difv = VIS*( (VO[ix-1][iy][iz]-2.0*VO[ix][iy][iz]+VO[ix+1][iy][iz])/(DX*DX)
                     +(VO[ix][iy-1][iz]-2.0*VO[ix][iy][iz]+VO[ix][iy+1][iz])/(DY*DY)
                     +(VO[ix][iy][iz-1]-2.0*VO[ix][iy][iz]+VO[ix][iy][iz+1])/(DZ*DZ) );
         /*���̑��x(V)�̌v�Z*/
         VN[ix][iy][iz] = VO[ix][iy][iz]
              + DT*( -cnvvx-cnvvy-cnvvz+difv+buov+(PO[ix][iy][iz]-PO[ix][iy+1][iz])/DY );
         }
       }
     }
     /*W(ix,iy,iz)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
         /* uu��W(ix,iy,iz)�ɂ�����U�̕�Ԓl */
         uu=(UO[ix-1][iy][iz]+UO[ix][iy][iz]+UO[ix-1][iy][iz+1]+UO[ix][iy][iz+1])/4.0;
         /* vv��W(ix,iy,iz)�ɂ�����V�̕�Ԓl */
         vv=(VO[ix][iy-1][iz]+VO[ix][iy][iz]+VO[ix][iy-1][iz+1]+VO[ix][iy][iz+1])/4.0;
         /* �Η���(cnvwx,cnvwy,cnvwz)���P�����x���㍷���ɂČv�Z */
         if ( uu >= 0.0 ){
           cnvwx = uu*( WO[ix][iy][iz] - WO[ix-1][iy][iz] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvwx = uu*( WO[ix+1][iy][iz] - WO[ix][iy][iz] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvwy = vv*( WO[ix][iy][iz]-WO[ix][iy-1][iz] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvwy = vv*( WO[ix][iy+1][iz]-WO[ix][iy][iz] ) / DY;
         }
         if ( WO[ix][iy][iz] >= 0.0 ){
           cnvwz = WO[ix][iy][iz]*( WO[ix][iy][iz]-WO[ix][iy][iz-1] ) / DZ;
         }
         else if ( WO[ix][iy][iz] < 0.0 ){
           cnvwz = WO[ix][iy][iz]*( WO[ix][iy][iz+1]-WO[ix][iy][iz] ) / DZ;
         }
         /*���͍�(buow)�̌v�Z*/
         tw = 0.0;
         buow = BUO*tw;
         /*�g�U��(difw)�̌v�Z*/
         difw = VIS*( (WO[ix-1][iy][iz]-2.0*WO[ix][iy][iz]+WO[ix+1][iy][iz])/(DX*DX)
                     +(WO[ix][iy-1][iz]-2.0*WO[ix][iy][iz]+WO[ix][iy+1][iz])/(DY*DY)
                     +(WO[ix][iy][iz-1]-2.0*WO[ix][iy][iz]+WO[ix][iy][iz+1])/(DZ*DZ) );
         /*���̑��x(W)�̌v�Z*/
         WN[ix][iy][iz] = WO[ix][iy][iz]
              + DT*( -cnvwx-cnvwy-cnvwz+difw+buow+(PO[ix][iy][iz]-PO[ix][iy][iz+1])/DZ );
         }
       }
     }
     velbnd();
}
/*���͏�̌v�Z Fortran �v���O�����̃R�����g���Q�Ƃ��Ă��������D*/
void press()
{
     /* --- SMAC --- ���͕␳�̐��`�V�X�e���p */
     double a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1],a6[NE1],a7[NE1],b[NE1],x[NE1];
     int ix,iy,iz,k,nnx,nny,nnz,nne,nnxy;
     /* ���`�V�X�e���֌W�̔z��̏����� --- SMAC --- */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         a1[k]=0.0; a2[k]=0.0; a3[k]=0.0; a4[k]=0.0; a5[k]=0.0; a6[k]=0.0; a7[k]=0.0;
         b[k]=0.0; x[k]=0.0;
         }
       }
     }
     /*P(ix,iy,iz)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           /* --- SMAC ---*/
           DIV[ix][iy][iz] = ( UN[ix][iy][iz] - UN[ix-1][iy][iz] )/DX
                           + ( VN[ix][iy][iz] - VN[ix][iy-1][iz] )/DY
                           + ( WN[ix][iy][iz] - WN[ix][iy][iz-1] )/DZ;
         }
       }
     }
     /* �W���s��̍쐬 --- SMAC ---
        �v�Z�̈��1�_����ɓ����̓_�i���z�Z�����܂߂čl�����2�_�����������̓_�j
        �Ɋւ��ČW���s����쐬�D�c��͋��E�����𔽉f������ pdbndc �Őݒ肷�� */
     for (ix=2;ix<=NX-1;ix++){
       for (iy=2;iy<=NY-1;iy++){
         for (iz=2;iz<=NZ-1;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[k] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ));
         a4[k] = - 1.0 / (DY*DY) * DT;
         a2[k] = - 1.0 / (DY*DY) * DT;
         a1[k] = - 1.0 / (DX*DX) * DT;
         a5[k] = - 1.0 / (DX*DX) * DT;
         a6[k] = - 1.0 / (DZ*DZ) * DT;
         a7[k] = - 1.0 / (DZ*DZ) * DT;
         b[k] = -DIV[ix][iy][iz];
         x[k] = PD[ix][iy][iz];
         }
       }
     }
     /* --- SMAC --- ���E�������W���s��ɔ��f������ */
     pdbndc(a1,a2,a3,a4,a5,a6,a7,b);
     /* --- SMAC --- ���͕␳ P'(PD) �Ɋւ���|�A�\���������̉�@ */
     nnx = NX; nny = NY; nnz= NZ; nne = NE; nnxy = NXY;
     /* 1. ���ږ@ : �o���h�}�g���b�N�X�ɂ��K�E�X�̏����@ */
     if (METHOD==1) gb    (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 2. �����@1 : point-SOR �@ */
     if (METHOD==2) psorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 3. �����@2 : line-SOR �@ : OMG=1�ŏ\���D�傫����������Ɣ��U���� */
     if (METHOD==3) lsorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 4. �N�����t������Ԗ@1 : �����c���@ */
     if (METHOD==4) crb   (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* 5. �N�����t������Ԗ@2 : BiCGSTAB */
     if (METHOD==5) bicgb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     /* METHOD=4,5�ɂ�����T���x�N�g���v�Z���̃[�����Z�΍� */
     if (IFLG==2)   psorb (a1,a2,a3,a4,a5,a6,a7,b,x,nnx,nny,nnz,nne,nnxy);
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         k = iy + (ix-1)*NY + (iz-1)*NXY;
         /* ���͂̑��ΐ��̏��� : ��l�̐ݒ�
            1���Ɨ��ȉ������߂�ꍇ�͈��͂̊�_�������I�ɐݒ肳��� */
         if (IRELP != 2) PD[ix][iy][iz] = x[k]; /* ���͂̊�_��݂��Ȃ��ꍇ */
         if (IRELP == 2) PD[ix][iy][iz] = x[k]-x[1]; /* P'(1,1,1)=0 ---> P(1,1,1)=0 */
         }
       }
     }
     /* ���͕␳�̋��E�����̏��� */
     pdbnd();
     /* ���x�̏C�� */
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
          UN[ix][iy][iz]=UN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix+1][iy][iz])/DX*DT;
         }
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         for (iz=1;iz<=NZ;iz++){
          VN[ix][iy][iz]=VN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix][iy+1][iz])/DY*DT;
         }
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ-1;iz++){
          WN[ix][iy][iz]=WN[ix][iy][iz]+(PD[ix][iy][iz]-PD[ix][iy][iz+1])/DZ*DT;
         }
       }
     }
     /*�V���ɓ���ꂽ���x��p���ċ��E��������������*/
     velbnd();
     /* ���͂̏C�� */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
          PN[ix][iy][iz] = PO[ix][iy][iz] + PD[ix][iy][iz];
         }
       }
     }
     /* �V���ɓ���ꂽ���͂�p���ċ��E�������������� */
     pbnd();
}
/*���x��̌v�Z*/
void caltem()
{
     int ix,iy,iz;
     double uut,vvt,wwt,cnvtx,cnvty,cnvtz,dift;
     /*T(ix,iy,iz)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
         /* uut,vvt,wwt�͂��ꂼT(ix,iy,iz)�ɂ�����U,V,W�̕�Ԓl */
         uut = ( UO[ix][iy][iz] + UO[ix-1][iy  ][iz  ] ) / 2.0;
         vvt = ( VO[ix][iy][iz] + VO[ix  ][iy-1][iz  ] ) / 2.0;
         wwt = ( WO[ix][iy][iz] + WO[ix  ][iy  ][iz-1] ) / 2.0;
         /*�Η���(cnvtx,cnvty,cnvtz)���P�����x���㍷���ɂČv�Z*/
         if ( uut >= 0.0 ){
           cnvtx = uut*( TO[ix][iy][iz] - TO[ix-1][iy][iz] ) / DX;
         }
         else if ( uut < 0.0 ){
           cnvtx = uut*( TO[ix+1][iy][iz] - TO[ix][iy][iz] ) / DX;
         }
         if ( vvt >= 0.0 ){
           cnvty = vvt*( TO[ix][iy][iz] - TO[ix][iy-1][iz] ) / DY;
         }
         else if ( vvt < 0.0 ){
           cnvty = vvt*( TO[ix][iy+1][iz] - TO[ix][iy][iz] ) / DY;
         }
         if ( wwt >= 0.0 ){
           cnvtz = wwt*( TO[ix][iy][iz] - TO[ix][iy][iz-1] ) / DZ;
         }
         else if ( wwt < 0.0 ){
           cnvtz = wwt*( TO[ix][iy][iz+1] - TO[ix][iy][iz] ) / DZ;
         }
         /* �g�U��(dift)�̌v�Z */
         dift = ALP*( (TO[ix-1][iy][iz]-2.0*TO[ix][iy][iz]+TO[ix+1][iy][iz])/(DX*DX)
                     +(TO[ix][iy-1][iz]-2.0*TO[ix][iy][iz]+TO[ix][iy+1][iz])/(DY*DY)
                     +(TO[ix][iy][iz-1]-2.0*TO[ix][iy][iz]+TO[ix][iy][iz+1])/(DZ*DZ) );
         /* ���̎��Ԃ�T�̌v�Z */
         TN[ix][iy][iz] = TO[ix][iy][iz] + DT*( -cnvtx-cnvty-cnvtz+dift );
         }
       }
     }
     /* ���E�����̏��� */
     tbnd();
}
/* ���x�̋��E�����̏��� */
void velbnd()
{
     int ix,iy,iz;
     for (iy=1;iy<=NY;iy++){/*U�i�E���ʁj*/
       for (iz=1;iz<=NZ;iz++){
         UN[NX][iy][iz] = 0.0;
       }
     }
     for (iy=1;iy<=NY;iy++){/*U�i�����ʁj*/
       for (iz=1;iz<=NZ;iz++){
         UN[0][iy][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX;ix++){/*U�i��ʁj*/
       for (iy=1;iy<=NY;iy++){
         UN[ix][iy][NZ+1] = -UN[ix][iy][NZ];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U�i�O�ʁj*/
       for (iy=1;iy<=NY;iy++){
         UN[ix][iy][0] = -UN[ix][iy][1];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U�i��ʁj*/
       for (iz=0;iz<=NZ+1;iz++){
         UN[ix][NY+1][iz] = -UN[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*U�i���ʁj*/
       for (iz=0;iz<=NZ+1;iz++){
         UN[ix][0][iz] = -UN[ix][1][iz];
       }
     }
     for (iy=1;iy<=NY-1;iy++){/*V�i�E���ʁj*/
       for (iz=1;iz<=NZ;iz++){
         VN[NX+1][iy][iz] = -VN[NX][iy][iz];
       }
     }
     for (iy=1;iy<=NY-1;iy++){/*V�i�����ʁj*/
       for (iz=1;iz<=NZ;iz++){
         VN[0][iy][iz] = -VN[1][iy][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i��ʁj*/
       for (iz=0;iz<=NZ+1;iz++){
         VN[ix][NY][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i���ʁj*/
       for (iz=0;iz<=NZ+1;iz++){
         VN[ix][0][iz] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i��ʁj*/
       for (iy=1;iy<=NY-1;iy++){
         VN[ix][iy][NZ+1] = -VN[ix][iy][NZ];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i�O�ʁj*/
       for (iy=1;iy<=NY-1;iy++){
         VN[ix][iy][0] = -VN[ix][iy][1];
       }
     }
     for (iy=1;iy<=NY;iy++){/*W�i�E�ʁj*/
       for (iz=1;iz<=NZ-1;iz++){
          WN[NX+1][iy][iz] = -WN[NX][iy][iz];
       }
     }
     for (iy=1;iy<=NY;iy++){/*W�i���ʁj*/
       for (iz=1;iz<=NZ-1;iz++){
          WN[0][iy][iz] = -WN[1][iy][iz];
       }
     }
     for (ix=0;ix<=NX;ix++){/*W�i��ʁj*/
       for (iy=0;iy<=NY+1;iy++){
          WN[ix][iy][NZ] = 0.0;
       }
     }
     for (ix=0;ix<=NX;ix++){/*W�i�O�ʁj*/
       for (iy=0;iy<=NY+1;iy++){
         WN[ix][iy][0] = 0.0;
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W�i��ʁj*/
       for (iz=1;iz<=NZ-1;iz++){
         WN[ix][NY+1][iz] = -WN[ix][NY][iz];
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*W�i���ʁj*/
       for (iz=1;iz<=NZ-1;iz++){
         WN[ix][0][iz] = -WN[ix][1][iz];
       }
     }
}
/*���x�̋��E�����̏���*/
void tbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*�E����*/
        for (iz=0;iz<=NZ+1;iz++){
          TN[NX+1][iy][iz] = 2.0 * ( -0.5 ) - TN[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*������*/
        for (iz=0;iz<=NZ+1;iz++){
          TN[0][iy][iz] = 2.0 * ( +0.5 ) - TN[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*���*/
        for (iz=1;iz<=NZ;iz++){
          TN[ix][NY+1][iz] = TN[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*����*/
        for (iz=1;iz<=NZ;iz++){
          TN[ix][0][iz] = TN[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*���*/
        for (iy=0;iy<=NY+1;iy++){
          TN[ix][iy][NZ+1] = TN[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/*�O��*/
        for (iy=0;iy<=NY+1;iy++){
          TN[ix][iy][0] = TN[ix][iy][1];
        }
      }
}
/* ���͕␳�̌W���s��Ƌ��E����
   �{�v�Z�v���O�����ł́C���͕␳�̋��E�������C�W���s��쐬���ɍl������D
   �܂��C����1���Ɨ����̊m�ۂ������ŏ�������D*/
void pdbndc(double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double a6[NE1], double a7[NE1], double b[NE1] )
{
     int ix,iy,iz,m;
     /* --- SMAC --- �W���s��ɋ��E�����𔽉f������
     �v�Z�̈�̍����̌W���s��i���E�����𔽉f�j*/
     ix = 1;
     for (iy=2;iy<=NY-1;iy++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̉E���̌W���s��i���E�����𔽉f�j*/
     ix = NX;
     for (iy=2;iy<=NY-1;iy++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) +1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̉����̌W���s��i���E�����𔽉f�j*/
     iy = 1;
     for (ix=2;ix<=NX-1;ix++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̏㑤�̌W���s��i���E�����𔽉f�j*/
     iy = NY;
     for (ix=2;ix<=NX-1;ix++){
       iz=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) +1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iz=2;iz<=NZ-1;iz++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iz=NZ;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̌�ʂ̌W���s��i���E�����𔽉f�j*/
     iz = 1;
     for (ix=2;ix<=NX-1;ix++){
       iy=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iy=2;iy<=NY-1;iy++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iy=NY;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a7[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̑O�ʂ̌W���s��i���E�����𔽉f�j*/
     iz = NZ;
     for (ix=2;ix<=NX-1;ix++){
       iy=1;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       for (iy=2;iy<=NY-1;iy++){
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
         a4[m] = - 1.0 / (DY*DY) * DT;
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
       }
       iy=NY;
         m = iy + (ix-1)*NY + (iz-1)*NXY;
         a3[m] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
         a2[m] = - 1.0 / (DY*DY) * DT;
         a1[m] = - 1.0 / (DX*DX) * DT;
         a5[m] = - 1.0 / (DX*DX) * DT;
         a6[m] = - 1.0 / (DZ*DZ) * DT;
         b[m] = -DIV[ix][iy][iz];
     }
     /* �v�Z�̈�̍���_��̌W���s��i���E�����𔽉f�j*/
     ix=1; iy=NY;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* �v�Z�̈�̍����_��̌W���s��i���E�����𔽉f�j*/
     ix=1; iy=1;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a5[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* �v�Z�̈�̉E��_��̌W���s��i���E�����𔽉f�j*/
     ix=NX; iy=NY;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a2[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* �v�Z�̈�̉E���_��̌W���s��i���E�����𔽉f�j*/
     ix=NX; iy=1;
     iz=1;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     for (iz=2;iz<=NZ-1;iz++){
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 2.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       a7[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     }
     iz=NZ;
       m = iy + (ix-1)*NY + (iz-1)*NXY;
       a3[m] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[m] = - 1.0 / (DY*DY) * DT;
       a1[m] = - 1.0 / (DX*DX) * DT;
       a6[m] = - 1.0 / (DZ*DZ) * DT;
       b[m] = -DIV[ix][iy][iz];
     /* �P���Ɨ��ȉ��𓾂邽�߂̏��� : ���ږ@�ɂ����Ă͕K�{
        (ix=1,iy=1,iz ---> m=1����_�Ƃ��C���PN[1][1][1]=PD[1][1][1]=0�Ƃ���) */
     if (IRELP ==1) {
       a3[1] = 1.0; a4[1] = 0.0; a5[1] = 0.0; a7[1]=0.0; b[1] = 0.0;
       /* m=2 �̓_�̏���(m=1�Ƃ̃����N��f��) */
       a3[2] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[2] = - 1.0 / (DY*DY) * DT;
       a2[2] = 0.0;
       a5[2] = - 1.0 / (DX*DX) * DT;
       a7[2] = - 1.0 / (DZ*DZ) * DT;
       /* m=1+NY �̓_�̏���(m=1�Ƃ̃����N��f��) */
       a3[1+NY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[1+NY] = - 1.0 / (DY*DY) * DT;
       a1[1+NY] = 0.0;
       a5[1+NY] = - 1.0 / (DX*DX) * DT;
       a7[1+NY] = - 1.0 / (DZ*DZ) * DT;
       /* m=1+NX*NY �̓_�̏���(m=1�Ƃ̃����N��f��) */
       a3[1+NXY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) + 1.0/(DZ*DZ) );
       a4[1+NXY] = - 1.0 / (DY*DY) * DT;
       a5[1+NXY] = - 1.0 / (DX*DX) * DT;
       a6[1+NXY] = 0.0;
       a7[1+NXY] = - 1.0 / (DZ*DZ) * DT;
     }
}
/********************************************************************
*                                                                   *
*                  ���͕␳�̋��E�����̏���                         *
*                                                                   *
* �{�v�Z�v���O�����ł́C�W���s��쐬���ɋ��E�������l�����Ă���̂ŁC*
* �����I�ɂ́C�ŏI�I�ɉ��z�Z���̒l��z��Ɋi�[���Ă���ɂ����Ȃ��D  *
* ���̃��[�`���͂Ȃ��Ƃ��悢�D                                      *
* �����ߒ��ŋ��E�������l������ꍇ�͕K�{�ŁC�����̂��тɂ��̏�����  *
* �s���K�v��������D                                                *
*                                                                   *
********************************************************************/
void pdbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*�E����*/
        for (iz=0;iz<=NZ+1;iz++){
          PD[NX+1][iy][iz] = PD[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*������*/
        for (iz=1;iz<=NZ+1;iz++){
          PD[0][iy][iz] = PD[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*���*/
        for (iz=0;iz<=NZ;iz++){
          PD[ix][NY+1][iz] = PD[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*����*/
        for (iz=1;iz<=NZ;iz++){
          PD[ix][0][iz] = PD[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/* ��� */
        for (iy=1;iy<=NY;iy++){
          PD[ix][iy][NZ+1] = PD[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/* �O�� */
        for (iy=1;iy<=NY;iy++){
          PD[ix][iy][0] = PD[ix][iy][1];
        }
      }
}
/*********************************************************************
*                     ���͂̋��E�����̏���
*********************************************************************/
void pbnd()
{
     int ix,iy,iz;
      for (iy=0;iy<=NY+1;iy++){/*�E����*/
        for (iz=0;iz<=NZ+1;iz++){
          PN[NX+1][iy][iz] = PN[NX][iy][iz];
        }
      }
      for (iy=0;iy<=NY+1;iy++){/*������*/
        for (iz=1;iz<=NZ+1;iz++){
          PN[0][iy][iz] = PN[1][iy][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*���*/
        for (iz=0;iz<=NZ;iz++){
          PN[ix][NY+1][iz] = PN[ix][NY][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/*����*/
        for (iz=1;iz<=NZ;iz++){
          PN[ix][0][iz] = PN[ix][1][iz];
        }
      }
      for (ix=1;ix<=NX;ix++){/* ��� */
        for (iy=1;iy<=NY;iy++){
          PN[ix][iy][NZ+1] = PN[ix][iy][NZ];
        }
      }
      for (ix=1;ix<=NX;ix++){/* �O�� */
        for (iy=1;iy<=NY;iy++){
          PN[ix][iy][0] = PN[ix][iy][1];
        }
      }
}
void prout()/*�f�[�^�o�͗p*/
{
    fwrite(UN, sizeof(double), NX1*NY2*NZ2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1*NZ2, out_12);
    fwrite(WN, sizeof(double), NX2*NY2*NZ1, out_13);
    fwrite(PN, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(PD, sizeof(double), NX2*NY2*NZ2, out_14);
    fwrite(TN, sizeof(double), NX2*NY2*NZ2, out_15);
}
/* Tecplot�p�f�[�^�̏o�� */
void tecplt ( char *fname_tec )
{

     int ix,iy,iz;
     double x,y,z,u,v,w,t;

     out_21=fopen(fname_tec,"wt");
     fprintf (out_21," VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"T\"\n");
     fprintf (out_21," ZONE I=%3d, J=%3d, K=%3d, F=POINT \n",NX1,NY1,NZ1);

     for (iz=0;iz<=NZ;iz++){
       for (iy=0;iy<=NY;iy++){
         for (ix=0;ix<=NX;ix++){
           x = DX * (double)ix;
           y = DY * (double)iy;
           z = DZ * (double)iz;
           u = (UN[ix][iy][iz]+UN[ix][iy+1][iz]+UN[ix][iy][iz+1]+UN[ix][iy+1][iz+1])/4.0;
           v = (VN[ix][iy][iz]+VN[ix+1][iy][iz]+VN[ix][iy][iz+1]+VN[ix+1][iy][iz+1])/4.0;
           w = (WN[ix][iy][iz]+WN[ix+1][iy][iz]+WN[ix][iy+1][iz]+WN[ix+1][iy+1][iz])/4.0;
           t = (TN[ix][iy][iz]+TN[ix+1][iy][iz]+TN[ix][iy+1][iz]+TN[ix+1][iy+1][iz]
               +TN[ix][iy][iz+1]+TN[ix+1][iy][iz+1]+TN[ix][iy+1][iz+1]+TN[ix+1][iy+1][iz+1] )/8.0;
         fprintf(out_21,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",x,y,z,u,v,w,t);
         }
       }
     }
     fclose (out_21);
}
/********************************************************************
*                    �e��̐��`�V�X�e����@                         *
*                                                                   *
*  1. �K�E�X�̏����@                                                *
*  2. point-SOR �@                                                  *
*  3. line-SOR �@                                                   *
*  4. �����c���@                                                    *
*  5. Bi-CGSTAB�@                                                   *
*                                                                   *
*  ��������C3�����̃|�A�\����������7�_�����ߎ��ɂė��U������       *
*  ���`�V�X�e�����������߂̂��̂ŁC�œK�����Ă���                   *
*  ������̃T�u���[�`������������Ƃ��Ă���                         *
*                                                                   *
********************************************************************/
/********************************************************************
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
********************************************************************/
void gb (double a1[NE1], double a2[NE1], double a3[NE1],
         double a4[NE1], double a5[NE1], double a6[NE1], double a7[NE1],
         double b[NE1], double x[NE1],
         int nx, int ny, int nz, int ne, int nxy)
{
     double a[NXYNXY][NE1];
     int ine,i,j,n,k;
     double aa,s;
     /* ���ږ@�̂Ƃ���(�W���s�񂪓��قłȂ����)�����Ȃ��ŁC�K�����𓾂� */
     IFLG = 0;
     /* �}�g���b�N�XA�̃[���N���A */
     for (ine=1;ine<=ne;ine++){
       for (i=-nxy;i<=nxy;i++){
         a[i+nxy][ine] = 0.0;
       }
     }
     /* �K�v�ȂƂ����A1����A7�܂ł��i�[���� */
     for (ine=1;ine<=ne;ine++){
       a[ -ny+nxy][ine] = a1[ine];
       a[  -1+nxy][ine] = a2[ine];
       a[   0+nxy][ine] = a3[ine];
       a[   1+nxy][ine] = a4[ine];
       a[  ny+nxy][ine] = a5[ine];
       a[-nxy+nxy][ine] = a6[ine];
       a[ nxy+nxy][ine] = a7[ine];
     }
/* �O�i���� */
     for (i=1;i<=ne-1;i++){
       if ( i <= ne-nxy ){
         for (j=1;j<=nxy;j++){
           aa = a[-j+nxy][i+j]/a[0+nxy][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=nxy-j;k++){
             a[k+nxy][i+j] = a[k+nxy][i+j]-a[n+nxy][i]*aa;
             n=n+1;
           }
         }
       }
       else{
         for (j=1;j<=ne-i;j++){
           aa = a[-j+nxy][i+j]/a[0+nxy][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ne-i-j;k++){
             a[k+nxy][i+j] = a[k+nxy][i+j]-a[n+nxy][i]*aa;
             n=n+1;
           }
         }
       }
     }
/* ��ޑ�� */
     /* �W���s��̓��ِ��𔻒� */
     if ( fabs( a[0+nxy][ne] ) <= 1.0e-50 ){
        printf("  Matrix is singular : |A(0,NE)| < 1E-50 \n");
        IFLG = 1; /* ���͕␳�̌v�Z�𔭎U�ɐݒ� */
     }
     x[ne] = b[ne] / a[0+nxy][ne];
     for (i=ne-1;i>=1;i--){
       s = 0.0;
       if ( i >  ne-nxy ){
         for (n=1;n<=ne-i;n++){
           s = s + a[n+nxy][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+nxy][i];
       }
       else{
         for (n=1;n<=nxy;n++){
           s = s + a[n+nxy][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+nxy][i];
       }
     }
}
/********************************************************************
*                                                                   *
*  point-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*  3�������v���V�A�����U���ɂ��7�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
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
********************************************************************/
void psorb (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double a6[NE1], double a7[NE1],
            double b[NE1], double x[NE1],
            int nx, int ny, int nz, int ne, int nxy)
{
     int i,j;
     double rnorm1, xold, sum, xnew;

     IFLG=1;/* FLG�̏����l��"��������" */

     for (j=1;j<=NITR;j++){
       rnorm1 = 0.0;

       i=1;
         xold = x[i];
         sum =                           +a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );

       for (i=2;i<=ny;i++){
         xold = x[i];
         sum =               a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ny+1;i<=nxy;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]               +a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=nxy+1;i<=ne-nxy;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a6[i]*x[i-nxy]+a7[i]*x[i+nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ne-nxy+1;i<=ne-ny;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny]+a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       for (i=ne-ny+1;i<=ne-1;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]              +a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }

       i=ne;
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]                           +a6[i]*x[i-nxy];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
 
       if ( rnorm1 <= EPSP ){/* ��������; �����Ȃ� IFLG=0 �ɐݒ� */
         IFLG=0; ITR = j; goto label_700;
       }
     }
 
     label_700:{};/* �����Ɣ��肳�ꂽ�Ƃ��̕���_ */
}
/********************************************************************
*                                                                   *
*  line-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`�� *
*  3�������v���V�A�����U���ɂ��7�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> at1[NE1],at2[NE1],at3[NE1],at4[NE1],at5[NE1],     *
*                 at6[NE1],at7[NE1]                                 *
*                                                                   *
*    B -> bx[NE1]: ���m�x�N�g��                                     *
*    X -> xn[NE1]: ���m�x�N�g�� ---> ��������߂�                   *
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
********************************************************************/
void lsorb (double at1[NE1] ,double at2[NE1], double at3[NE1],
            double at4[NE1], double at5[NE1], double at6[NE1], double at7[NE1],
            double bx[NE1], double xn[NE1],
            int nx, int ny, int nz, int ne, int nxy)
{
     double x1[NE1],xo[NE1],xold[NE1],a[NE1],b[NE1],c[NE1],d[NE1],u[NE1],y[NE1];
     double rnorm;
     int k, inx, iy, ix, iz, iny, j, i, inz, inzz;

     IFLG=1; /* IFLG �̏����l��"��������" */

     for (k=1;k<=NITR;k++){
       inx = 1;/* x �������ւ̑|�� : �g�[�}�X�@�ɂ��*/
       for (iz=1;iz<=nz;iz++){
         for (iy=1;iy<=ny;iy++){
           for (ix=1;ix<=nx;ix++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X */
           a[inx]=at1[inz]; b[inx]=at3[inz]; c[inx]=at5[inz]; d[inx]=bx[inz];
           if (inz-1 >= 1){
             d[inx]=d[inx]-at2[inz]*xn[inz-1];
           }
           if (inz+1 <= ne){
             d[inx]=d[inx]-at4[inz]*xn[inz+1];
           }
           if (inz-nxy >= 1){
             d[inx]=d[inx]-at6[inz]*xn[inz-nxy];
           }
           if (inz+nxy <= ne){
             d[inx]=d[inx]-at7[inz]*xn[inz+nxy];
           }
           xold[inz]=xn[inz];/* x�����ւ̑|���̑O�̒l��XOLD�ɕۑ� */
           xo[inz]=xn[inz];/* �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ� */
           inx=inx+1;
           }
         }
       }
       u[1]=c[1]/b[1];/* Ly=b ������ */
       for (j=2;j<=ne-1;j++){
         u[j]=c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1]=d[1]/b[1];
       for (j=2;j<=ne;j++){
         y[j]=( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y ������ */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       inx = 1;
       for (iz=1;iz<=nz;iz++){
         for (iy=1;iy<=ny;iy++){
           for (ix=1;ix<=nx;ix++){
           inz=iy+(ix-1)*ny+(iz-1)*nxy;
           /* ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa */
           xn[inz]=(1.0-OMG)*xo[inz]+OMG*x1[inx];
           inx=inx+1;
           }
         }
       }
       iny = 1;/* y �������ւ̑|�� : �g�[�}�X�@�ɂ�� */
       for (iz=1;iz<=nz;iz++){
         for (ix=1;ix<=nx;ix++){
           for (iy=1;iy<=ny;iy++){
           inz = iy + (ix-1)*nx + (iz-1)*nxy;
           /* �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X */
           a[iny]=at2[inz]; b[iny]=at3[inz]; c[iny]=at4[inz]; d[iny]=bx[inz];
           if (inz-ny >= 1){
             d[iny]=d[iny]-at1[inz]*xn[inz-ny];
           }
           if (inz+ny <= ne){
             d[iny]=d[iny]-at5[inz]*xn[inz+ny];
           }
           if (inz-nxy >= 1){
             d[iny]=d[iny]-at6[inz]*xn[inz-nxy];
           }
           if (inz+nxy <= ne){
             d[iny]=d[iny]-at7[inz]*xn[inz+nxy];
           }
           xo[inz]=xn[inz];/* �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ� */
           iny=iny+1;
           }
         }
       }
       u[1] = c[1] / b[1];/* Ly=b ������ */
       for (j=2;j<=ne-1;j++){
         u[j] = c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         y[j] = ( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y ������ */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       iny = 1;
       for (iz=1;iz<=nz;iz++){
         for (ix=1;ix<=nx;ix++){
           for (iy=1;iy<=ny;iy++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa */
           xn[inz] = (1.0-OMG)*xo[inz]+OMG*x1[iny];
           iny = iny + 1;
           }
         }
       }
       inzz = 1;/* z �������ւ̑|�� : �g�[�}�X�@�ɂ�� */
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           for (iz=1;iz<=nz;iz++){
           inz = iy + (ix-1)*nx + (iz-1)*nxy;
           /* �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X */
           a[inzz]=at6[inz]; b[inzz]=at3[inz]; c[inzz]=at7[inz]; d[inzz]=bx[inz];
           if (inz-ny >= 1){
             d[inzz]=d[inzz]-at1[inz]*xn[inz-ny];
           }
           if (inz+ny <= ne){
             d[inzz]=d[inzz]-at5[inz]*xn[inz+ny];
           }
           if (inz-1 >= 1){
             d[inzz]=d[inzz]-at2[inz]*xn[inz-1];
           }
           if (inz+1 <= ne){
             d[inzz]=d[inzz]-at4[inz]*xn[inz+1];
           }
           xo[inz]=xn[inz];/* �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ� */
           inzz=inzz+1;
           }
         }
       }
       u[1] = c[1] / b[1];/* Ly=b ������ */
       for (j=2;j<=ne-1;j++){
         u[j] = c[j]/( b[j]-a[j]*u[j-1] );
       }
       y[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         y[j] = ( d[j]-a[j]*y[j-1] ) / ( b[j]-a[j]*u[j-1] );
       }
       x1[ne] = y[ne];/* Ux=y ������ */
       for (j=ne-1;j>=1;j--){
         x1[j] = y[j] - u[j]*x1[j+1];
       }
       inzz = 1;
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           for (iz=1;iz<=nz;iz++){
           inz = iy + (ix-1)*ny + (iz-1)*nxy;
           /* ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa */
           xn[inz] = (1.0-OMG)*xo[inz]+OMG*x1[inzz];
           inzz = inzz + 1;
           }
         }
       }
       rnorm= 0.0;
       for (i=1;i<=ne;i++){
         rnorm= rnorm + (xn[i]-xold[i])*(xn[i]-xold[i]);
       }
       if (rnorm <= EPSP){/* �������� */
         IFLG=0; ITR=k; goto label_900;
       }
     }
     label_900:{};/* �����Ɣ��肳�ꂽ�Ƃ��̕���_ */
}
/********************************************************************
*                                                                   *
*    �����c��(Conjugate Residual)�@�ɂ���Ώ̍s�� A ���܂�        *
*  ���`�V�X�e����@�T�u���[�`��                                     *
*                        AX=B                                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
*                                                                   *
* B : ���m�x�N�g��                                                  *
* X : ���m�x�N�g�� ---> ��������߂� ---> �����ł͕֋X��z�� XP[NE1]*
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
*    R[NE1]: r_{k} = B - A x_{k}                                    *
*    P[NE1]:  p_{k+1} = r_{k+1} + ��_{k} p_{k}, p_{0} = r_{0}       *
*    AP[NE1]: A * P                                                 *
*    AR[NE1]: A * R                                                 *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
********************************************************************/
void crb (double a1[NE1], double a2[NE1], double a3[NE1],
          double a4[NE1], double a5[NE1], double a6[NE1], double a7[NE1],
          double b[NE1], double xp[NE1],
          int nx, int ny, int nz, int ne, int nxy)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1], x[NE1], xold[NE1];/* ��Ɨp�z�� */
     double rap, apap, alp, rnorm, beta, arap;
     int i, k;

     promv (a1,a2,a3,a4,a5,a6,a7,x,r,nx,ny,nz,ne,nxy);/* R �� AX ���� */
     for (i=1;i<=ne;i++){
       r[i]=b[i]-r[i]; p[i]=r[i];/* r_{0}��p_{0}(�����l)�̐ݒ� */
       xold[i]=xp[i];/* �O�̎�����X��XOLD�ɑ�� */
     }
     promv (a1,a2,a3,a4,a5,a6,a7,p,ap,nx,ny,nz,ne,nxy);/* AP�� A p_{0} ���� */
/* �����v�Z */
     for (k=1;k<=NITR;k++){
       rap = provv (r,ap,ne);/* ( r_{k}, A p_{k} )�̌v�Z => RAP */
       apap = provv (ap,ap,ne);/* ( A p_{k}, A p_{k} )�̌v�Z => APAP */
       if (fabs(apap) < 1.0e-50){
         printf(" 0 division : ALPHA_{K} in Conjugate Residual \n");
         IFLG = 2; return;/* ���`�V�X�e���̌v�Z�𔭎U�ɐݒ� */
       }
       else {
         alp = rap / apap;
       }
       rnorm = 0.0;
       for (i=1;i<=ne;i++){
         x[i]=x[i]+alp*p[i];/* x_{k+1}=x_{k}+��_{k}p_{k} */
         r[i]=r[i]-alp*ap[i];/* r_{k+1}=r_{k}-��_{k}Ap_{k} */
         /* �O�̔����Ƃ̍��̃m�����̌v�Z */
         rnorm=rnorm+(x[i]-xold[i])*(x[i]-xold[i]);
         xold[i]=x[i];/* ����ꂽX��XOLD�ɑ�� */
       }
       /* rnorm �� EPSP �ȉ��Ȃ�����Ƃ݂Ȃ��� 700 �� */
       if (rnorm <= EPSP){
         IFLG=0; ITR=k; goto label_700;
       }
       /* ���������̏ꍇ */
       promv (a1,a2,a3,a4,a5,a6,a7,r,ar,nx,ny,nz,ne,nxy);/* A r_{k+1} �̌v�Z => AR(NE) */
       arap =  provv( ar,ap,ne);/* ( A r_{k+1}, A p_{k} )�̌v�Z => ARAP */
       if (fabs(apap) < 1.0e-50){
         printf(" 0 division : BETA_{K} in Conjugate Residual \n");
         IFLG = 2; return;/* ���`�V�X�e���̌v�Z�𔭎U�ɐݒ� */
       }
       else {
         beta = -arap / apap;/* ��_{k} = - ( A r_{k+1}, A p_{k} ) / ( A p_{k}, A p_{k} ) */
       }
       for (i=1;i<=ne;i++){
         p[i]=r[i]+beta*p[i];/* p_{k+1} = r_{k+1} + ��_{k} p_{k} */
         ap[i]=ar[i]+beta*ap[i];/* A p_{k+1} = A r_{k+1} + ��_{k}A p_{k} */
       }
     }
     IFLG=1;/* NITR �܂Ōv�Z���Ă��������� */
     label_700:{}; /* �����Ɣ��肳�ꂽ�Ƃ��̕���_ */
     for (i=1;i<=ne;i++){
       xp[i]=x[i];
     }
}
/********************************************************************
*                                                                   *
*    �x�N�g�� A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��                 *
*                        AB=C                                       *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NE : ���i�q�_��(�x�N�g�� A,B �̃T�C�Y)                         *
*    C  : A �� B �̐�(�X�J���[)                                     *
*                                                                   *
********************************************************************/
double provv (double a[NE1], double b[NE1], int ne)
{
      int i;
      double c;

      c = 0.0;
      for (i=1;i<=ne;i++){
        c = c + a[i]*b[i];
      }
      return(c);
}
/********************************************************************
*                                                                   *
*    �}�g���b�N�X A �ƃx�N�g�� B �̐ς̌v�Z�T�u���[�`��             *
*                        AB=C                                       *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NE : ���i�q�_��(�����}�g���b�N�X A,B,C �̃T�C�Y)               *
*    C  : A �� B �̐�(�x�N�g��)                                     *
*                                                                   *
********************************************************************/
void promv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
            double a5[NE1], double a6[NE1], double a7[NE1],
            double b[NE1], double c[NE1],
            int nx, int ny, int nz, int ne, int nxy)
{
     int i;

     i=1; c[i]=a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];

     for (i=2;i<=ny;i++){
       c[i]=a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];
     }

     for (i=ny+1;i<=nxy;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]+a7[i]*b[i+nxy];
     }

     for (i=nxy+1;i<=ne-nxy;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]
            +a6[i]*b[i-nxy]+a7[i]*b[i+nxy];
     }

     for (i=ne-nxy+1;i<=ne-ny;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny]
           +a6[i]*b[i-nxy];
     }
     for (i=ne-ny+1;i<=ne-1;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a6[i]*b[i-nxy];
     }

     i=ne; c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a6[i]*b[i-nxy];

}
/********************************************************************
*                                                                   *
*  Bi-CGSTAB �@�ɂ���Ώ̍s�� A ���܂�                            *
*  ���`�V�X�e����@�T�u���[�`��                                     *
*                        AX=B                                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],                  *
*                 a5[NE1],a6[NE1],a7[NE1]                           *
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
* [�z��̐���]                                                      *
*    T[NE] : t_{k} = r_{k} - ��_{k} A p_{k}                         *
*    X[NE] : x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{K} t_{k}          *
*    R[NE] : r_{k+1} = t_{k} - ��_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P[NE] : p_{k+1} = r_{k+1} + ��_{k} ( p_{k}-��_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP[NE] : A * P                                                 *
*    AR[NE] : A * T                                                 *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
********************************************************************/
void bicgb ( double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
             double a6[NE1], double a7[NE1], double b[NE1], double x[NE1],
             int nx, int ny, int nz, int ne, int nxy)
{
     double r[NE1], ap[NE1], at[NE1], p[NE1], s[NE1], t[NE1], xold[NE1];
     int i,j;
     double sr1, sap, alpha, att, atat, xi, sr2, beta, rnorm;

     for (i=1;i<=ne;i++){
       xold[i] = x[i];
     }

      promv (a1,a2,a3,a4,a5,a6,a7,x,r,nx,ny,nz,ne,nxy);/* R �� AX ���� */
      for (i=1;i<=ne;i++){/* r_{0}��p_{0}(�����l)�C������ s=r_{0} �̐ݒ� */
        r[i] = b[i] - r[i]; p[i] = r[i]; s[i] = r[i];
      }

      /* �J��Ԃ��v�Z */
      for (j=1;j<=NITR;j++){
        sr1=provv (s,r,ne);/* ( s, r_{k} ) �̌v�Z => SR1 */
        promv (a1,a2,a3,a4,a5,a6,a7,p,ap,nx,ny,nz,ne,nxy);/* A p_{k} �̌v�Z => AP(NE) */
        sap=provv (s,ap,ne);/* ( s, A p_{k} ) �̌v�Z => SAP */
        if (fabs(sap) < 1.0e-50){
          printf(" 0 division : ALPHA_{K} in Bi-CGSTAB \n");
          IFLG = 2; return;/* ���`�V�X�e���̌v�Z�𔭎U�ɐݒ� */
        }
        else {
          alpha = sr1/sap;/* ��_{k} = ( s, r_{k} ) / ( s, A p_{k} ) */
        }
        for (i=1;i<=ne;i++){
          t[i] = r[i] - alpha*ap[i];/* t_{k} = r_{k} - ��_{k} A p_{k} */
        }
        promv (a1,a2,a3,a4,a5,a6,a7,t,at,nx,ny,nz,ne,nxy);/* A t_{k} �̌v�Z => AT(NE) */
        att=provv (at,t,ne);/* ( A t_{k}, t_{k} ) �̌v�Z => ATT */
        atat=provv (at,at,ne);/* ( A t_{k}, A t_{k} ) �̌v�Z => ATAT */
        if (fabs(atat) < 1.0e-50){
          printf(" 0 division : XI_{K} in Bi-CGSTAB \n");
          IFLG = 2; return;/* ���`�V�X�e���̌v�Z�𔭎U�ɐݒ� */
        }
        else {
          xi = att / atat;/* ��_{k} = ( A t_{k}, t_{k} ) / ( A t_{k}, A t_{k} ) */
        }
        rnorm = 0.0;
        for (i=1;i<=ne;i++){
          x[i] = x[i] + alpha*p[i] + xi*t[i];/* x_{k+1}=x_{k}+��_{k}p_{k}+��_{k}t_{k} */
          r[i] = t[i] - xi*at[i];/* r_{k+1}=t_{k}-��_{k} A t_{k} */
          /* �O�̔����Ƃ̍��̃m�����̌v�Z */
          rnorm = rnorm + (x[i]-xold[i])*(x[i]-xold[i]);
          xold[i] = x[i];/* ����ꂽX��XOLD�ɑ�� */
        }
        if (rnorm <= EPSP){/* RNORM �� EPSP �ȉ��Ȃ�����Ƃ݂Ȃ��� 900 �� */
          IFLG=0; ITR=j; goto label_900;
        }
        /* ���������̏ꍇ */
        sr2=provv (s,r,ne);/* ( s, r_{k+1} ) �̌v�Z => SR2 */
        /* ��_{k} = ( ��_{k} / ��_{k} ) * ( s, r_{k+1} )/( s, r_{k} ) */
        beta = (alpha / xi) * (sr2 / sr1);
        for (i=1;i<=ne;i++){
          /* p_{k+1} = r_{k+1} + ��_{k}( p_{k}-��_{k} A p_{k} ) */
          p[i] = r[i] + beta * ( p[i] - xi*ap[i] );
        }
      }
      /* NITR �܂Ōv�Z���Ă��������� */
      IFLG=1;

      label_900:{};/* �����Ɣ��肳�ꂽ�Ƃ��̕���_ */
}
