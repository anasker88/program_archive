/********************************************************************
* �t�@�C����  �Fhsmac3d.c                                           *
* �^�C�g��    �FHSMAC�@�ɂ��3�����M������̓v���O����              *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� ���p���w��                      *
* �����      �F2003.12.25                                          *
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
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NZ1 21 /* NZ1=NZ+1 */
#define NZ2 22 /* NZ2=NZ+2 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
double DX, DY, DZ, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
char fname[12][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_15, *out_21;
FILE *in_16, *in_17, *in_18, *in_19, *in_20;
double UO[NX1][NY2][NZ2],UN[NX1][NY2][NZ2],VO[NX2][NY1][NZ2],VN[NX2][NY1][NZ2];
double WO[NX2][NY2][NZ1],WN[NX2][NY2][NZ1],PO[NX2][NY2][NZ2],TO[NX2][NY2][NZ2],TN[NX2][NY2][NZ2];
void main()
{
     int i,j,NITR,NCYCLE;
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
     printf(" DLX = %12.3lE  DLY = %12.3lE  DLZ = %12.3lE  IRELP = %d  METHOD = %d \n",
              DLX,DLY,DLZ,IRELP,METHOD);
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
     calvel();/*���x��̌v�Z*/
     ITR=1;
     label_710:{};/*���͔����̂��߂̖߂�_*/
     IFLG=0;
     press();/*���̈��͏�̌v�Z*/
     if ( IFLG==0 ){/*Newton-Raphson�@�ɂ�鈳�͏�̌v�Z�����������Ƃ�*/
       if ( ITYPE==2 ){/*�񓙉���v�Z�̏ꍇ*/
         caltem();/*���x����v�Z*/
       }
     }
     if ( IFLG==1 ){/*���͏�̌v�Z���������Ă��Ȃ��Ƃ�*/
       if ( ITR<NITR ){/*���͌v�Z�̔����񐔂����炩���ߐݒ肳�ꂽ�ő�lNITR��菬�����Ƃ�*/
         ITR = ITR + 1;/*����ɔ������J��Ԃ�*/
         goto label_710;
       }
       else{/*���͌v�Z�̔����񐔂�NITR�ȏ�̂Ƃ����U�Ƃ݂Ȃ��Čv�Z�I��*/
         printf(" calculation has diverged \n");
         prout();
         goto label_900;
       }
     }
     if ( ICYCLE < NCYCLE ){/*���Ԑi�s�T�C�N��(ICYCLE)��NCYCLE��菬������*/
       goto label_700;
     }
     else{/*���Ԑi�s�T�C�N����NCYCLE�ɂȂ�����->�v�Z�I��*/
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
             PO[ix][iy][iz] = 0.0;
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
       fread(PO, sizeof(double), NX2*NY2*NZ2, in_19);
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
         for (iz=0;iz<=NZ+1;iz++){/*UO:�V�������ԃX�e�b�v�ł̏����l�DUN��ۑ��D*/
           UO[ix][iy][iz] = UN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*    VN -> VO : �K�v�Ȃ����ւ���O��VN��VO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY;iy++){/*    VN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         for (iz=0;iz<=NZ+1;iz++){/*VO:�V�������ԃX�e�b�v�ł̏����l�DVN��ۑ��D*/
           VO[ix][iy][iz] = VN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/*WN -> WO : �K�v�Ȃ����ւ���O��WN��WO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*WN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         for (iz=0;iz<=NZ;iz++){/*WO:�V�������ԃX�e�b�v�ł̏����l�DWN��ۑ��D*/
           WO[ix][iy][iz] = WN[ix][iy][iz];
         }
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : �K�v�Ȃ����ւ���O��TN��TO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*  TN:�O�̎��ԃX�e�b�v�ł̌v�Z�l*/
         for (iz=0;iz<=NZ+1;iz++){/*TO:�V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D*/
           TO[ix][iy][iz] = TN[ix][iy][iz];
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
/*���͏�̌v�Z*/
void press()
{
     int ixmax,iymax,izmax,ix,iy,iz;
     double del,div,delp,postn;
     ixmax = 0; iymax = 0; izmax=0; DMAX = 0.0e0;
     /*P(ix,iy,iz)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         for (iz=1;iz<=NZ;iz++){
           del = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) + 2.0/(DZ*DZ) );
           div = ( UN[ix][iy][iz] - UN[ix-1][iy][iz] )/DX
               + ( VN[ix][iy][iz] - VN[ix][iy-1][iz] )/DY
               + ( WN[ix][iy][iz] - WN[ix][iy][iz-1] )/DZ;
           if ( fabs(div) >= fabs(DMAX) ){
             ixmax = ix; iymax = iy; izmax = iz; DMAX = div;
           }
           delp = - OMG * div / del;
           PO[ix][iy][iz] = PO[ix][iy][iz] + delp;
           UN[ix  ][iy  ][iz  ]=UN[ix  ][iy  ][iz  ]+DT/DX*delp;
           UN[ix-1][iy  ][iz  ]=UN[ix-1][iy  ][iz  ]-DT/DX*delp;
           VN[ix  ][iy  ][iz  ]=VN[ix  ][iy  ][iz  ]+DT/DY*delp;
           VN[ix  ][iy-1][iz  ]=VN[ix  ][iy-1][iz  ]-DT/DY*delp;
           WN[ix  ][iy  ][iz  ]=WN[ix  ][iy  ][iz  ]+DT/DZ*delp;
           WN[ix  ][iy  ][iz-1]=WN[ix  ][iy  ][iz-1]-DT/DZ*delp;
         }
       }
     }
     /* ���͂̑��ΐ��Ɋւ��鏈��(IRELP=1�Ȃ�ȉ��̏������s��) */
     if (IRELP==1){
       postn = PO[1][1][1];
       for (ix=1;ix<=NX;ix++){
         for (iy=1;iy<=NY;iy++){
           for (iz=1;iz<=NZ;iz++){
              PO[ix][iy][iz] = PO[ix][iy][iz]-postn;
           }
         }
       }
     }
     /* IFLG=1�Ȃ�C�A���̎��𖞂����Ă��Ȃ��Ɣ��肵�Ăш��͌v�Z���s�� */
     if ( fabs(DMAX) >= EPSP ){
       IFLG = 1;
     }
     /* ���͌v�Z�̉񐔂�100�񂲂Ƃɕ\�� */
     if ( (ITR%100) == 0){
       printf (" Iteration= %8d, Div(max)( %6d, %6d, %6d ) = %15.6lE\n",ITR,ixmax,iymax,izmax,DMAX);
     }
     /*�V���ɓ���ꂽ���x��p���ċ��E��������������*/
     velbnd();
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
void prout()/*�f�[�^�o�͗p*/
{
    fwrite(UN, sizeof(double), NX1*NY2*NZ2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1*NZ2, out_12);
    fwrite(WN, sizeof(double), NX2*NY2*NZ1, out_13);
    fwrite(PO, sizeof(double), NX2*NY2*NZ2, out_14);
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
