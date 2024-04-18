/********************************************************************
* �t�@�C����  �Fsmac2d.c                                            *
* �^�C�g��    �FSMAC�@�ɂ��2�����M������̓v���O����               *
* �����      �F����@���V                                          *
* ����        �F���R���ȑ�w �H�w�� �o�C�I�E���p���w��              *
* �����      �F2011.11.01                                          *
* ����        �FC                                                   *
*********************************************************************
* �v���O���������s����ƁC"in2d.mac"(�K���������j��ǂ݂ɍs���D     *
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
#define NE 400 /* NE:���͕␳�̖��m���̐� */
#define NX1 21 /* NX1=NX+1 */
#define NX2 22 /* NX2=NX+2 */
#define NY1 21 /* NY1=NY+1 */
#define NY2 22 /* NY2=NY+2 */
#define NYNY 41 /* NYNY=2*NY+1 */
#define NE1 401 /* NE1=NE+1 */
void cinit();
void adv();
void calvel();
void press();
void caltem();
void velbnd();
void tbnd();
void prout();
void tecplt( char *fname_tec);
void pdbndc (double [NE1], double [NE1], double [NE1], double [NE1], double [NE1], double [NE1]);
void gb (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
         double [NE1], double [NE1], int nx, int ny, int ne);
void crb   (double [NE1], double [NE1], double [NE1], double[NE1], double [NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
void promv (double [NE1], double [NE1], double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1],
            int nx, int ny, int ne);
void bicgb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
double provv (double [NE1], double [NE1], int ne);
void psorb (double [NE1], double [NE1], double [NE1], double[NE1], double[NE1],
            double [NE1], double [NE1], int nx, int ny, int ne);
void lsorb (double [NE1] ,double [NE1], double [NE1],
            double [NE1], double [NE1], double [NE1], double [NE1],
           int nx, int ny, int ne);
void pdbnd();
void pbnd();
double DX, DY, DT;
double VIS, ALP, BUO;
double RE, PR, GR, TIME, OMG, EPSP;
int  ICYCLE, ITR, IFLG, IRELP, METHOD;
double DMAX;
int ITYPE;
int NITR;
char fname[10][80];
char *fname10;
FILE *in_10, *out_11, *out_12, *out_13, *out_14, *out_21;
FILE *in_15, *in_16, *in_17, *in_18;
double UO[NX1][NY2], UN[NX1][NY2], VO[NX2][NY1], VN[NX2][NY1];
double PO[NX2][NY2], PN[NX2][NY2], TO[NX2][NY2], TN[NX2][NY2];
double DIV[NX1][NY1], PD[NX2][NY2];/* --- SMAC ---*/
void main()
{
     int i,j,NCYCLE;
     double DLX,DLY;
     char buff[80];
     fname10="in2d.mac";
     strncpy(fname[0],fname10,8);
     in_10=fopen(fname[0],"rt"); /*�p�����[�^�t�@�C���̃I�[�v��*/
     for (i=1;i<=9;i++){/*�o�̓t�@�C�����̓ǂݍ���*/
       fgets(buff, sizeof buff,in_10);
       for (j=0;buff[j]!=0;j++);
       strncpy(fname[i], buff, j-1);
     }
     out_11=fopen(fname[1],"wb");/*U�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_12=fopen(fname[2],"wb");/*V�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_13=fopen(fname[3],"wb");/*P�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     out_14=fopen(fname[4],"wb");/*T�̌v�Z���ʏo�͗p�t�@�C���I�[�v��*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(10�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*                   ��(11�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(12�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(13�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(14�s��)*/
     fgets(buff, sizeof buff,in_10);/*                   ��(15�s��)*/
     fscanf(in_10," %d %d %d %d ",&ITYPE,&ICYCLE,&NITR,&NCYCLE);
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(17�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(18�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg ",&EPSP,&OMG);
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(20�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(21�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg %lg %lg ",&DT,&RE,&PR,&GR);
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(23�s��)�̃X�L�b�v*/
     fgets(buff, sizeof buff,in_10);/*in2d.mac���̃R�����g�s(24�s��)�̃X�L�b�v*/
     fscanf(in_10," %lg %lg %d %d ",&DLX,&DLY,&IRELP,&METHOD);
     printf (" ITYPE = %d ICYCLE= %d NITR= %d NCYCLE= %d \n",
             ITYPE,ICYCLE,NITR,NCYCLE);
     printf(" EPSP= %12.3lE  OMG = %12.3lE \n",EPSP,OMG);
     printf(" DT = %12.3lE  RE = %12.3lE  PR = %12.3lE  GR = %12.3lE \n",DT,RE,PR,GR);
     printf(" DLX = %12.3lE  DLY = %12.3lE IRELP = %d  METHOD = %d \n",DLX,DLY,IRELP,METHOD);
     if ( ICYCLE != 0 ){/*�p���̌v�Z�̏ꍇ*/
       in_15=fopen(fname[5],"rb");/*U�f�[�^�t�@�C���̃I�[�v��*/
       in_16=fopen(fname[6],"rb");/*V�f�[�^�t�@�C���̃I�[�v��*/
       in_17=fopen(fname[7],"rb");/*P�f�[�^�t�@�C���̃I�[�v��*/
       in_18=fopen(fname[8],"rb");/*T�f�[�^�t�@�C���̃I�[�v��*/
     }
     DX=DLX/(double)NX; DY=DLY/(double)NY;/*x,y�����̊i�q��DX,DY*/
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
     else{/*���Ԑi�s�J�E���^��NCYCLE�ɂȂ�����v�Z�I��*/
       prout();
     }
     label_900:{};
     tecplt(fname[9]);/*Tecplot�p�f�[�^�̏o��*/
     fclose (in_10); fclose (out_11); fclose (out_12);
     fclose (out_13); fclose (out_14);
}
/*�����ݒ�*/
void cinit()
{
     int ix,iy;
     if ( ICYCLE == 0 ){/*�V�K�v�Z�̏ꍇ*/
       for (ix=0;ix<=NX;ix++){/*U�̏����l�ݒ�*/
         for (iy=0;iy<=NY+1;iy++){
           UN[ix][iy] = 0.0;
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*V�̏����l�ݒ�*/
         for (iy=0;iy<=NY;iy++){
           VN[ix][iy] = 0.0;
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*P�̏����l�ݒ�*/
         for (iy=0;iy<=NY+1;iy++){
           PD[ix][iy] = 0.0; PN[ix][iy] = 0.0;/* --- SMAC ---*/
         }
       }
       for (ix=0;ix<=NX+1;ix++){/*T�̏����l�ݒ�(�̈���͍���(+0.5)�ƒቷ(-0.5)�̒��ԉ��x)*/
         for (iy=0;iy<=NY+1;iy++){/*�i���Ӂj���͍��̌v�Z�ŉ��x�̔z����g�p����*/
           TN[ix][iy] = 0.0;/*����̂œ�����ł�T=0�Ƃ��ď������������͐ݒ肷��*/
         }/*�K�v������D�[���ȊO�̒l������ƕ��͍����v�Z�����\��������̂Œ��ӁD*/
       }
       for (iy=0;iy<=NY+1;iy++){/*T�̋��E�F�E���ǁi��p�jT=-0.5*/
         TN[NX+1][iy] = 2.0 * ( -0.5 ) - TN[NX][iy];
       }
       for (iy=0;iy<=NY+1;iy++){/*T�̋��E�F�����ǁi���M�jT=+0.5*/
         TN[0][iy] = 2.0 * ( +0.5 ) - TN[1][iy];
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F��ʁi�f�M�j*/
         TN[ix][NY+1] = TN[ix][NY];
       }
       for (ix=1;ix<=NX;ix++){/*T�̋��E�F���ʁi�f�M�j*/
         TN[ix][0] = TN[ix][1];
       }
     }
     else{/*�p���v�Z�i���łɂ���v�Z���ʂ���X�^�[�g�j�̏ꍇ*/
       /*U�f�[�^�t�@�C������̓ǂݍ���[Unit No.=15]*/
       fread(UN, sizeof(double), NX1*NY2, in_15);
       /*V�f�[�^�t�@�C������̓ǂݍ���[Unit No.=16]*/
       fread(VN, sizeof(double), NX2*NY1, in_16);
       /*P�f�[�^�t�@�C������̓ǂݍ���[Unit No.=17]*/
       fread(PN, sizeof(double), NX2*NY2, in_17);
       fread(PD, sizeof(double), NX2*NY2, in_17);
       /*T�f�[�^�t�@�C������̓ǂݍ���[Unit No.=18]*/
       fread(TN, sizeof(double), NX2*NY2, in_18);
       fclose (in_15); fclose (in_16); fclose (in_17); fclose (in_18);
     }
}
/*���Ԑi�s*/
void adv()
{
     int ix,iy;
     TIME = DT*(double)ICYCLE; ICYCLE = ICYCLE + 1;
     if ( (ICYCLE%100) ==0 ){/*ICYCLE��100�񖈂ɕ\��*/
       printf ("  CYC = %8d \n",ICYCLE);
     }
     for (ix=0;ix<=NX;ix++){/* UN -> UO : �K�v�Ȃ����ւ���O��UN��UO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*UN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         UO[ix][iy] = UN[ix][iy];/*UO:�V�������ԃX�e�b�v�ł̏����l�DUN��ۑ��D*/
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* VN -> VO : �K�v�Ȃ����ւ���O��VN��VO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY;iy++){/*VN:�O�̎��ԃX�e�b�v�ɂ����čŏI�I�ɓ���ꂽ�l�C���͕␳�̓x�ɍX�V�����*/
         VO[ix][iy] = VN[ix][iy];/*VO:�V�������ԃX�e�b�v�ł̏����l�DVN��ۑ��D*/
       }
     }
     for (ix=0;ix<=NX+1;ix++){/* TN -> TO : �K�v�Ȃ����ւ���O��TN��TO����ϓ��ʂ����߂�*/
       for (iy=0;iy<=NY+1;iy++){/*TN:�O�̎��ԃX�e�b�v�ł̌v�Z�l*/
         TO[ix][iy] = TN[ix][iy];/*TO:�V�������ԃX�e�b�v�ł̏����l�DTN��ۑ��D*/
         PO[ix][iy] = PN[ix][iy];/* --- SMAC --- PO:�V�������ԃX�e�b�v�ł̏����l�DPN��ۑ��D*/
       }
     }
}
/*���x��̌v�Z*/
void calvel()
{
     int ix,iy;
     double vv,cnvux,cnvuy,tu,buou,difu;
     double uu,cnvvx,cnvvy,tv,buov,difv;
     /*U(ix,iy)�̌v�Z*/
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
         /* vv��U(ix,iy)�ɂ�����V�̕�Ԓl */
         vv=(VO[ix][iy]+VO[ix+1][iy]+VO[ix][iy-1]+VO[ix+1][iy-1])/4.0;
         /* �Η���(cnvux,cnvuy)���P�����x���㍷���ɂČv�Z */
         if ( UO[ix][iy] >= 0.0 ){
           cnvux = UO[ix][iy]*( UO[ix][iy]-UO[ix-1][iy] ) / DX;
         }
         else if ( UO[ix][iy] < 0.0 ){
           cnvux = UO[ix][iy]*( UO[ix+1][iy] - UO[ix][iy] ) / DX;
         }
         if ( vv >= 0.0 ){
           cnvuy = vv*( UO[ix][iy] - UO[ix][iy-1] ) / DY;
         }
         else if ( vv < 0.0 ){
           cnvuy = vv*( UO[ix][iy+1] - UO[ix][iy] ) / DY;
         }
         /* x�����̕��͍�(buou)�̓[�� */
         tu = 0.0;
         buou = 0.0;
         /* �g�U��(difu)�̌v�Z */
         difu = VIS*( (UO[ix-1][iy]-2.0*UO[ix][iy]+UO[ix+1][iy])/(DX*DX)
                     +(UO[ix][iy-1]-2.0*UO[ix][iy]+UO[ix][iy+1])/(DY*DY) );
         /*���̑��x(U)�̌v�Z*/
         UN[ix][iy] = UO[ix][iy]
              + DT*( -cnvux-cnvuy+difu+buou+( PO[ix][iy]-PO[ix+1][iy] )/DX );
       }
     }
     /*V(ix,iy)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
         /* uu��V(ix,iy)�ɂ�����U�̕�Ԓl */
         uu=(UO[ix-1][iy]+UO[ix][iy]+UO[ix-1][iy+1]+UO[ix][iy+1])/4.0;
         /* �Η���(cnvvx,cnvvy)���P�����x���㍷���ɂČv�Z */
         if ( uu >= 0.0 ){
           cnvvx = uu*( VO[ix][iy] - VO[ix-1][iy] ) / DX;
         }
         else if ( uu < 0.0 ){
           cnvvx = uu*( VO[ix+1][iy] - VO[ix][iy] ) / DX;
         }
         if ( VO[ix][iy] >= 0.0 ){
           cnvvy = VO[ix][iy]*( VO[ix][iy]-VO[ix][iy-1] ) / DY;
         }
         else if ( VO[ix][iy] < 0.0 ){
           cnvvy = VO[ix][iy]*( VO[ix][iy+1]-VO[ix][iy] ) / DY;
         }
         /*���͍�(buov)�̌v�Z*/
         tv = ( TO[ix][iy] + TO[ix][iy+1] )/2.0;
         buov = BUO*tv;
         /*�g�U��(difv)�̌v�Z*/
         difv = VIS*( (VO[ix-1][iy]-2.0*VO[ix][iy]+VO[ix+1][iy])/(DX*DX)
                     +(VO[ix][iy-1]-2.0*VO[ix][iy]+VO[ix][iy+1])/(DY*DY) );
         /*���̑��x(V)�̌v�Z*/
         VN[ix][iy] = VO[ix][iy]
              + DT*( -cnvvx-cnvvy+difv+buov+(PO[ix][iy]-PO[ix][iy+1])/DY );
       }
     }
     velbnd();/* ���x�̋��E�����̏��� */
}
/*���͏�̌v�Z Fortran �v���O�����̃R�����g���Q�Ƃ��Ă��������D*/
void press()
{
     double a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1],b[NE1],x[NE1];/* --- SMAC --- ���͕␳�̐��`�V�X�e���p */
     int ix,iy,k,nnx,nny,nne;
     /* ���`�V�X�e���֌W�̔z��̏����� :*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         k = iy + (ix-1)*NY;
         a1[k]=0.0; a2[k]=0.0; a3[k]=0.0; a4[k]=0.0; a5[k]=0.0;
         b[k]=0.0; x[k]=0.0;
       }
     }
     /*P(ix,iy)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         /* --- SMAC ---*/
         DIV[ix][iy] = ( UN[ix][iy] - UN[ix-1][iy  ] )/DX
                     + ( VN[ix][iy] - VN[ix  ][iy-1] )/DY;
       }
     }
     /* �W���s��̍쐬
        �v�Z�̈��1�_����ɓ����̓_�i���z�Z�����܂߂čl�����2�_�����������̓_�j
        �Ɋւ��ČW���s����쐬�D�c��͋��E�����𔽉f������ pdbndc �Őݒ肷�� */
     for (ix=2;ix<=NX-1;ix++){
       for (iy=2;iy<=NY-1;iy++){
         k = iy + (ix-1)*NY;
         a3[k] = DT*( 2.0/(DX*DX) + 2.0/(DY*DY) );
         a4[k] = - 1.0 / (DY*DY) * DT;
         a2[k] = - 1.0 / (DY*DY) * DT;
         a1[k] = - 1.0 / (DX*DX) * DT;
         a5[k] = - 1.0 / (DX*DX) * DT;
         b[k] = -DIV[ix][iy];
         x[k] = PD[ix][iy];
       }
     }
     /* --- SMAC --- ���E�������W���s��ɔ��f������ */
     pdbndc(a1,a2,a3,a4,a5,b);
     /* --- SMAC --- ���͕␳ P'(PD) �Ɋւ���|�A�\���������̉�@ */
     nnx = NX; nny = NY; nne = NE;
     /* 1. ���ږ@ : �o���h�}�g���b�N�X�ɂ��K�E�X�̏����@ */
     if (METHOD==1) gb    (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 2. �����@1 : point-SOR �@ */
     if (METHOD==2) psorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 3. �����@2 : line-SOR �@ : OMG=1�ŏ\���D�傫����������Ɣ��U���� */
     if (METHOD==3) lsorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 4. �N�����t������Ԗ@1 : �����c���@ */
     if (METHOD==4) crb   (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* 5. �N�����t������Ԗ@2 : BiCGSTAB */
     if (METHOD==5) bicgb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* METHOD=4,5�ɂ�����T���x�N�g���v�Z���̃[�����Z�΍� */
     if (IFLG==2)   psorb (a1,a2,a3,a4,a5,b,x,nnx,nny,nne);
     /* ���`�V�X�e����@�œ���ꂽ1�����z��̉���2�����z��ɒu�������� */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         k=iy+(ix-1)*NY;
         /* ���͂̑��ΐ��̏��� : ��l�̐ݒ�
            1���Ɨ��ȉ������߂�ꍇ�͈��͂̊�_�������I�ɐݒ肳��� */
         if (IRELP != 2)  PD[ix][iy] = x[k]; /* ���͂̊�_��݂��Ȃ��ꍇ */
         /* 1���]���ȉ��̂�����1�����߂���C���͊��݂���ꍇ (IRELP=2) */
         if (IRELP == 2) PD[ix][iy] = x[k]-x[1]; /* P'(1,1)=0 ---> P(1,1)=0 */
       }
     }
     /* ���͕␳�̋��E�����̏��� */
     pdbnd();
     /* ���x�̏C�� */
     for (ix=1;ix<=NX-1;ix++){
       for (iy=1;iy<=NY;iy++){
          UN[ix][iy] = UN[ix][iy] + ( PD[ix][iy]-PD[ix+1][iy] )/DX*DT;
       }
     }
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY-1;iy++){
          VN[ix][iy] = VN[ix][iy] + ( PD[ix][iy]-PD[ix][iy+1] )/DY*DT;
       }
     }
     /*�V���ɓ���ꂽ���x��p���ċ��E��������������*/
     velbnd();
     /* ���͂̏C�� */
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
          PN[ix][iy] = PO[ix][iy] + PD[ix][iy];
       }
     }
     /* �V���ɓ���ꂽ���͂�p���ċ��E�������������� */
     pbnd();
}
/*���x��̌v�Z*/
void caltem()
{
     int ix,iy;
     double uut,vvt,cnvtx,cnvty,dift;
     /*T(ix,iy)�̌v�Z*/
     for (ix=1;ix<=NX;ix++){
       for (iy=1;iy<=NY;iy++){
         /* uut,vvt�͂��ꂼT(IX,IY)�ɂ�����U,V�̕�Ԓl */
         uut = ( UO[ix][iy] + UO[ix-1][iy  ] ) / 2.0;
         vvt = ( VO[ix][iy] + VO[ix  ][iy-1] ) / 2.0;
         /*�Η���(cnvtx,cnvty)���P�����x���㍷���ɂČv�Z*/
         if ( uut >= 0.0 ){
           cnvtx = uut*( TO[ix][iy] - TO[ix-1][iy] ) / DX;
         }
         else if ( uut < 0.0 ){
           cnvtx = uut*( TO[ix+1][iy] - TO[ix][iy] ) / DX;
         }
         if ( vvt >= 0.0 ){
           cnvty = vvt*( TO[ix][iy] - TO[ix][iy-1] ) / DY;
         }
         else if ( vvt < 0.0 ){
           cnvty = vvt*( TO[ix][iy+1] - TO[ix][iy] ) / DY;
         }
         /* �g�U��(dift)�̌v�Z */
         dift = ALP*( (TO[ix-1][iy]-2.0*TO[ix][iy]+TO[ix+1][iy])/(DX*DX)
                     +(TO[ix][iy-1]-2.0*TO[ix][iy]+TO[ix][iy+1])/(DY*DY) );
         /* ���̎��Ԃ�T�̌v�Z */
         TN[ix][iy] = TO[ix][iy] + DT*( -cnvtx-cnvty+dift );
       }
     }
     /* ���E�����̏��� */
     tbnd();
}
/* ���x�̋��E�����̏��� */
void velbnd()
{
     int ix,iy;
     for (iy=1;iy<=NY;iy++){/*U�i�E���ʁj*/
       UN[NX][iy] = 0.0;
     }
     for (iy=1;iy<=NY;iy++){/*U�i�����ʁj*/
       UN[0][iy] = 0.0;
     }
     for (ix=0;ix<=NX;ix++){/*U�i��ʁj*/
       UN[ix][NY+1] = -UN[ix][NY];
     }
     for (ix=0;ix<=NX;ix++){/*U�i���ʁj*/
       UN[ix][0] = -UN[ix][1];
     }
     for (iy=1;iy<=NY-1;iy++){/*V�i�E���ʁj*/
       VN[NX+1][iy] = -VN[NX][iy];
     }
     for (iy=1;iy<=NY-1;iy++){/*V�i�����ʁj*/
       VN[0][iy] = -VN[1][iy];
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i��ʁj*/
       VN[ix][NY] = 0.0;
     }
     for (ix=0;ix<=NX+1;ix++){/*V�i���ʁj*/
       VN[ix][0] = 0.0;
     }
}
/*���x�̋��E�����̏���*/
void tbnd()
{
     int ix,iy;
      for (iy=0;iy<=NY+1;iy++){/*�E����*/
        TN[NX+1][iy] = 2.0 * ( -0.5 ) - TN[NX][iy];
      }
      for (iy=0;iy<=NY+1;iy++){/*������*/
        TN[0][iy] = 2.0 * ( +0.5 ) - TN[1][iy];
      }
      for (ix=1;ix<=NX;ix++){/*���*/
        TN[ix][NY+1] = TN[ix][NY];
      }
      for (ix=1;ix<=NX;ix++){/*����*/
        TN[ix][0] = TN[ix][1];
      }
}
/* ���͕␳�̌W���s��Ƌ��E����
   �{�v�Z�v���O�����ł́C���͕␳�̋��E�������C�W���s��쐬���ɍl������D
   �܂��C����1���Ɨ����̊m�ۂ������ŏ�������D*/
void pdbndc(double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
            double b[NE1] )
{
     int ix,iy,k;
     /* --- SMAC --- �W���s��ɋ��E�����𔽉f������
     �v�Z�̈�̍����̌W���s��i���E�����𔽉f�j*/
     ix = 1;
     for (iy=2;iy<=NY-1;iy++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a2[k] = - 1.0 / (DY*DY) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* �v�Z�̈�̉E���̌W���s��i���E�����𔽉f�j*/
     ix = NX;
     for (iy=2;iy<=NY-1;iy++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a2[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* �v�Z�̈�̉����̌W���s��i���E�����𔽉f�j*/
     iy = 1;
     for (ix=2;ix<=NX-1;ix++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a4[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* �v�Z�̈�̏㑤�̌W���s��i���E�����𔽉f�j*/
     iy = NY;
     for (ix=2;ix<=NX-1;ix++){
       k = iy + (ix-1)*NY;
       a3[k] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a2[k] = - 1.0 / (DY*DY) * DT;
       a1[k] = - 1.0 / (DX*DX) * DT;
       a5[k] = - 1.0 / (DX*DX) * DT;
       b[k] = -DIV[ix][iy];
     }
     /* �����_�̌W���s��i���E�����𔽉f�j*/
     ix = 1; iy = 1; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a4[k] = - 1.0 / (DY*DY) * DT;
     a5[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* ����_�̌W���s��i���E�����𔽉f�j*/
     ix = 1; iy = NY; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a2[k] = - 1.0 / (DY*DY) * DT;
     a5[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* �E���_�̌W���s��i���E�����𔽉f�j*/
     ix = NX; iy = 1; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a4[k] = - 1.0 / (DY*DY) * DT;
     a1[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* �E��_�̌W���s��i���E�����𔽉f�j*/
     ix = NX; iy = NY; k = iy + (ix-1)*NY;
     a3[k] = DT*( 1.0/(DX*DX) + 1.0/(DY*DY) );
     a2[k] = - 1.0 / (DY*DY) * DT;
     a1[k] = - 1.0 / (DX*DX) * DT;
     b[k] = -DIV[ix][iy];
     /* �P���Ɨ��ȉ��𓾂邽�߂̏��� : ���ږ@�ɂ����Ă͕K�{
        (ix=1,iy=1 ---> k=1����_�Ƃ��C���PN[1][1]=PD[1][1]=0�Ƃ���) */
     if (IRELP==1) {
       a3[1] = 1.0; a4[1] = 0.0; a5[1] = 0.0; b[1] = 0.0;
       /* k=2 �̓_�̏���(k=1�Ƃ̃����N��f��) */
       a3[2] = DT*( 1.0/(DX*DX) + 2.0/(DY*DY) );
       a4[2] = - 1.0 / (DY*DY) * DT;
       a2[2] = 0.0;
       a5[2] = - 1.0 / (DX*DX) * DT;
       /* k=1+NY �̓_�̏���(k=1�Ƃ̃����N��f��) */
       a3[1+NY] = DT*( 2.0/(DX*DX) + 1.0/(DY*DY) );
       a4[1+NY] = - 1.0 / (DY*DY) * DT;
       a1[1+NY] = 0.0;
       a5[1+NY] = - 1.0 / (DX*DX) * DT;
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
     int ix,iy;
      for (iy=1;iy<=NY;iy++){/*�E����*/
        PD[NX+1][iy] = PD[NX][iy];
      }
      for (iy=1;iy<=NY;iy++){/*������*/
        PD[0][iy] = PD[1][iy];
      }
      for (ix=0;ix<=NX+1;ix++){/*���*/
        PD[ix][NY+1] = PD[ix][NY];
      }
      for (ix=0;ix<=NX+1;ix++){/*����*/
        PD[ix][0] = PD[ix][1];
      }
}
/*********************************************************************
*                     ���͂̋��E�����̏���
*********************************************************************/
void pbnd()
{
     int ix,iy;
      for (iy=1;iy<=NY;iy++){/*�E����*/
        PN[NX+1][iy] = PN[NX][iy];
      }
      for (iy=1;iy<=NY;iy++){/*������*/
        PN[0][iy] = PN[1][iy];
      }
      for (ix=0;ix<=NX+1;ix++){/*���*/
        PN[ix][NY+1] = PN[ix][NY];
      }
      for (ix=0;ix<=NX+1;ix++){/*����*/
        PN[ix][0] = PN[ix][1];
      }
}
void prout()/*�f�[�^�o�͗p*/
{
    fwrite(UN, sizeof(double), NX1*NY2, out_11);
    fwrite(VN, sizeof(double), NX2*NY1, out_12);
    fwrite(PN, sizeof(double), NX2*NY2, out_13);
    fwrite(PD, sizeof(double), NX2*NY2, out_13);
    fwrite(TN, sizeof(double), NX2*NY2, out_14);
}
/* Tecplot�p�f�[�^�̏o�� */
void tecplt ( char *fname_tec )
{

     int ix,iy;
     double x,y,u,v,t;

     out_21=fopen(fname_tec,"wt");
     fprintf (out_21," VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"T\"\n");
     fprintf (out_21," ZONE I=%3d, J=%3d, F=POINT \n",NX1,NY1);

     for (iy=0;iy<=NY;iy++){
       for (ix=0;ix<=NX;ix++){
         x = DX * (double)ix;
         y = DY * (double)iy;
         u = ( UN[ix][iy]+UN[ix][iy+1])/2.0;
         v = ( VN[ix][iy]+VN[ix+1][iy])/2.0;
         t = ( TN[ix][iy]+TN[ix+1][iy]
              +TN[ix][iy+1]+TN[ix+1][iy+1] )/4.0;
         fprintf(out_21,"%12.3lE %12.3lE %12.3lE %12.3lE %12.3lE\n",x,y,u,v,t);
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
*  ��������C2�����̃|�A�\����������5�_�����ߎ��ɂė��U������       *
*  ���`�V�X�e�����������߂̂��̂ŁC�œK�����Ă���                   *
*  ������̃T�u���[�`������������Ƃ��Ă���                         *
*                                                                   *
********************************************************************/
/********************************************************************
*                                                                   *
*   2�������v���V�A�����U���ɂ��5�_�����ߎ�                        *
*   �ɂē���ꂽ�K���I��Ώ̍s��A���܂񂾐��`�V�X�e��               *
*                        AX=B                                       *
*   ��Gauss�̏����@��p���ĉ����T�u���[�`���D(���I��)             *
*   A�̓o���h�}�g���b�N�X                                           *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*            ---> a[-NY:NY][NE1] -> a[NYNY][NE1]�Ɋi�[���Ȃ���      *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*                                                                   *
********************************************************************/
void gb (double a1[NE1], double a2[NE1], double a3[NE1],
         double a4[NE1], double a5[NE1], double b[NE1],
         double x[NE1], int nx, int ny, int ne)
{
     double a[NYNY][NE1];
     int ine,i,j,n,k;
     double aa,s;

     /* ���ږ@�̂Ƃ���(�W���s�񂪓��قłȂ����)�����Ȃ��ŁC�K�����𓾂� */
     IFLG = 0;

     /* �}�g���b�N�XA�̃[���N���A */
     for (ine=1;ine<=ne;ine++){
       for (i=-ny;i<=ny;i++){
         a[i+ny][ine] = 0.0;
       }
     }

     /* �K�v�ȂƂ����AT1����AT5�܂ł��i�[���� */
     for (ine=1;ine<=ne;ine++){
       a[-ny+ny][ine] = a1[ine];
       a[ -1+ny][ine] = a2[ine];
       a[  0+ny][ine] = a3[ine];
       a[  1+ny][ine] = a4[ine];
       a[ ny+ny][ine] = a5[ine];
     }

/* �O�i���� */
     for (i=1;i<=ne-1;i++){
       if ( i <= ne-ny ){
         for (j=1;j<=ny;j++){
           aa = a[-j+ny][i+j]/a[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ny-j;k++){
             a[k+ny][i+j] = a[k+ny][i+j]-a[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
       else{
         for (j=1;j<=ne-i;j++){
           aa = a[-j+ny][i+j]/a[0+ny][i];
           b[i+j] = b[i+j] - b[i]*aa;
           n=1;
           for (k=-j+1;k<=ne-i-j;k++){
             a[k+ny][i+j] = a[k+ny][i+j]-a[n+ny][i]*aa;
             n=n+1;
           }
         }
       }
     }

/* ��ޑ�� */
     /* �W���s��̓��ِ��𔻒� */
     if ( fabs( a[0+ny][ne] ) <= 1.0e-50 ){
        printf("  Matrix Singular : |A(0,NE)| < 1E-50 \n");
        IFLG = 1; return;/* ���͕␳�̌v�Z�𔭎U�ɐݒ� */
     }
     x[ne] = b[ne] / a[0+ny][ne];
     for (i=ne-1;i>=1;i--){
       s = 0.0;
       if ( i >  ne-ny ){
         for (n=1;n<=ne-i;n++){
           s = s + a[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+ny][i];
       }
       else{
         for (n=1;n<=ny;n++){
           s = s + a[n+ny][i] * x[i+n];
         }
         x[i] = ( b[i]-s ) / a[0+ny][i];
       }
     }

}
/********************************************************************
*                                                                   *
*  point-SOR �@�ɂ���Ώ̍s�� A ���܂ސ��`�V�X�e����@�T�u���[�`��*
*  2�������v���V�A�����U���ɂ��5�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in2d.mac)�ɂĐݒ�                *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
********************************************************************/
void psorb (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], 
            double a5[NE1], double b[NE1], double x[NE1],
            int nx, int ny, int ne)
{
     int i,j;
     double rnorm1, xold, sum, xnew;

     IFLG=1;/* FLG�̏����l��"��������" */
 
     for (j=1;j<=NITR;j++){
       rnorm1 = 0.0;
 
       i=1;
         xold = x[i];
         sum =                           +a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
 
       for (i=2;i<=ny;i++){
         xold = x[i];
         sum =               a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       for (i=ny+1;i<=ne-ny;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1]+a5[i]*x[i+ny];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       for (i=ne-ny+1;i<=ne-1;i++){
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1]+a4[i]*x[i+1];
         xnew = ( b[i]-sum )/a3[i];
         x[i] = xold + OMG * ( xnew - xold );
         rnorm1 = rnorm1 + ( xnew - xold )*( xnew - xold );
       }
 
       i=ne;
         xold = x[i];
         sum = a1[i]*x[i-ny]+a2[i]*x[i-1];
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
*  2�������v���V�A�����U���ɂ��5�_�����ߎ��p                       *
*                                                                   *
*    ���`�V�X�e�� --->  A_{i,j} X_{i} = B_{i}                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B -> bx[NE1] : ���m�x�N�g��                                    *
*    X -> xn[NE1] : ���m�x�N�g�� ---> ��������߂�                  *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in2d.mac)�ɂĐݒ�                *
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
            double at4[NE1], double at5[NE1], double bx[NE1], double xn[NE1],
            int nx, int ny, int ne)
{
     double x1[NE1],xo[NE1],xold[NE1],a[NE1],b[NE1],c[NE1],d[NE1],p[NE1],q[NE1];
     double rnorm;
     int k, inx, iy, ix, iny, j, i ;

     IFLG=1; /* IFLG �̏����l��"��������" */

     for (k=1;k<=NITR;k++){/* x �������ւ̑|�� : �g�[�}�X�@�ɂ��*/
       inx = 1;
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           iny = iy + (ix-1)*ny;
           /* �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X */
           a[inx]=at1[iny]; b[inx]=at3[iny]; c[inx]=at5[iny]; d[inx]=bx[iny];
           if (iny-1 >= 1){
             d[inx]=d[inx]-at2[iny]*xn[iny-1];
           }
           if (iny+1 <= ne){
             d[inx]=d[inx]-at4[iny]*xn[iny+1];
           }
           xold[iny]=xn[iny];/* x�����ւ̃g�[�}�X�@�œ��������߂�O��XN��XOLD�ɕۑ� */
           xo[iny]=xn[iny];/* �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ� */
           inx=inx+1;
         }
       }
       p[1]=c[1]/b[1];/* Ly=b ������ */
       for (j=2;j<=ne-1;j++){
         p[j]=c[j]/( b[j]-a[j]*p[j-1] );
       }
       q[1]=d[1]/b[1];
       for (j=2;j<=ne;j++){
         q[j]=( d[j]-a[j]*q[j-1] ) / ( b[j]-a[j]*p[j-1] );
       }
       x1[ne] = q[ne];/* Ux=y ������ */
       for (j=ne-1;j>=1;j--){
         x1[j] = q[j] - p[j]*x1[j+1];
       }
       inx = 1;
       for (iy=1;iy<=ny;iy++){
         for (ix=1;ix<=nx;ix++){
           iny=iy+(ix-1)*ny;
           /* ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa */
           xn[iny]=(1.0-OMG)*xo[iny]+OMG*x1[inx];
           inx=inx+1;
         }
       }
       iny = 1;/* y �������ւ̑|�� : �g�[�}�X�@�ɂ�� */
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           iny = ix + (iy-1)*nx;
           /* �g�[�}�X�@�̂��߂̌W��A,B,C,D�̐ݒ� : XN�͍ŐV��X */
           a[iny]=at2[iny]; b[iny]=at3[iny]; c[iny]=at4[iny]; d[iny]=bx[iny];
           if (iny-ny >= 1){
             d[iny]=d[iny]-at1[iny]*xn[iny-ny];
           }
           if (iny+ny <= ne){
             d[iny]=d[iny]-at5[iny]*xn[iny+ny];
           }
           xo[iny]=xn[iny];/* �g�[�}�X�@�œ��������߂�O��XN��XO�ɕۑ� */
           iny=iny+1;
         }
       }
       p[1] = c[1] / b[1];/* Ly=b ������ */
       for (j=2;j<=ne-1;j++){
         p[j] = c[j]/( b[j]-a[j]*p[j-1] );
       }
       q[1] = d[1] / b[1];
       for (j=2;j<=ne;j++){
         q[j] = ( d[j]-a[j]*q[j-1] ) / ( b[j]-a[j]*p[j-1] );
       }
       x1[ne] = q[ne];/* Ux=y ������ */
       for (j=ne-1;j>=1;j--){
         x1[j] = q[j] - p[j]*x1[j+1];
       }
       iny = 1;
       for (ix=1;ix<=nx;ix++){
         for (iy=1;iy<=ny;iy++){
           /* ����ꂽX1�Ɣ����O��XO�ɂ��ŐV��XN���ɘa */
           xn[iny] = (1.0-OMG)*xo[iny]+OMG*x1[iny];
           iny = iny + 1;
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
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*  B : ���m�x�N�g��                                                 *
*  X : ���m�x�N�g�� ---> ��������߂� ---> �����ł͕֋X��z��xp[NE1]*
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in2d.mac)�ɂĐݒ�                *
*                                                                   *
*�m�z��̐����n                                                     *
*    r[NE] : r_{k} = B - A x_{k}                                    *
*    p[NE] :  p_{k+1} = r_{k+1} + ��_{k} p_{k}, p_{0} = r_{0}       *
*    ap[NE] : A * P                                                 *
*    ar[NE] : A * R                                                 *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
********************************************************************/
void crb (double a1[NE1], double a2[NE1], double a3[NE1],
          double a4[NE1], double a5[NE1], double b[NE1],
          double xp[NE1], int nx, int ny, int ne)
{
     double r[NE1], p[NE1], ap[NE1], ar[NE1], x[NE1], xold[NE1];/* ��Ɨp�z�� */
     double rap, apap, alp, rnorm, beta, arap;
     int i, k;

     promv (a1,a2,a3,a4,a5,x,r,nx,ny,ne);/* R �� AX ���� */
     for (i=1;i<=ne;i++){
       r[i]=b[i]-r[i]; p[i]=r[i];/* r_{0}��p_{0}(�����l)�̐ݒ� */
       xold[i]=xp[i];/* �O�̎�����X��XOLD�ɑ�� */
     }
     promv (a1,a2,a3,a4,a5,p,ap,nx,ny,ne);/* AP�� A p_{0} ���� */
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
       promv (a1,a2,a3,a4,a5,r,ar,nx,ny,ne);/* A r_{k+1} �̌v�Z => AR(NE) */
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
*    ne : ���i�q�_��(�x�N�g�� A,B �̃T�C�Y)                         *
*    c  : a �� b �̐�(�X�J���[)                                     *
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
*    ne : ���i�q�_��(�����}�g���b�N�X A, �x�N�g��B,C �̃T�C�Y)      *
*    c  : a �� b �̐�(�x�N�g��)                                     *
*                                                                   *
********************************************************************/
void promv (double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1],
            double a5[NE1], double b[NE1], double c[NE1],
            int nx, int ny, int ne)
{
     int i;

     i=1; c[i]=a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];

     for (i=2;i<=ny;i++){
       c[i]=a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];
     }

     for (i=ny+1;i<=ne-ny;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1]+a5[i]*b[i+ny];
     }

     for (i=ne-ny+1;i<=ne-1;i++){
       c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i]+a4[i]*b[i+1];
     }

     i=ne; c[i]=a1[i]*b[i-ny]+a2[i]*b[i-1]+a3[i]*b[i];

}
/********************************************************************
*                                                                   *
*  Bi-CGSTAB �@�ɂ���Ώ̍s�� A ���܂�                            *
*  ���`�V�X�e����@�T�u���[�`��                                     *
*                        AX=B                                       *
*                                                                   *
*    �W���s��̌v�Z�e�ʐߖ�F�ڍׂ̓T�u���[�`�� PRESS ���Q��        *
*    A_{i,j} ---> a1[NE1],a2[NE1],a3[NE1],a4[NE1],a5[NE1]           *
*                                                                   *
*    B : ���m�x�N�g��                                               *
*    X : ���m�x�N�g�� ---> ��������߂�                             *
*                                                                   *
*�m�ϐ��̐����n                                                     *
*    NX : x�����i�q������                                           *
*    NY : y�����i�q������                                           *
*    NE : ���i�q�_�� = NX * NY                                      *
*    NITR : ���e������(in2d.mac)�ɂĐݒ�                          *
*    EPSP : ������������ŗp����l(in2d.mac)�ɂĐݒ�                *
*                                                                   *
*�m�z��̐����n                                                     *
*    T[NE] : t_{k} = r_{k} - ��_{k} A p_{k}                         *
*    X[NE] : x_{k+1} = x_{k} + ��_{k} p_{k} + ��_{K} t_{k}          *
*    R[NE] : r_{k+1} = t_{k} - ��_{k} A t_{k}                       *
*            r_{0} = B - A x_{0}                                    *
*    P[NE] : p_{k+1} = r_{k+1} + ��_{k} ( p_{k}-��_{k} A p_{k} )    *
*            p_{0} = r_{0}                                          *
*    AP[NE] : A * P                                                 *
*    AT[NE] : A * T                                                 *
*                                                                   *
* [�����������]                                                    *
*  (\vec{x}^{new}-\vec{x}^{old})^{2} < EPSP                         *
*  \vec{x}^{old} : �O�̔����ɂ��l                                 *
*  \vec{x}^{new} : �V���������ɂ��l                               *
*                                                                   *
********************************************************************/
void bicgb ( double a1[NE1], double a2[NE1], double a3[NE1], double a4[NE1], double a5[NE1],
             double b[NE1], double x[NE1], int nx, int ny, int ne)
{
     double r[NE1], ap[NE1], at[NE1], p[NE1], s[NE1], t[NE1], xold[NE1];
     int i,j;
     double sr1, sap, alpha, att, atat, xi, sr2, beta, rnorm;

     for (i=1;i<=ne;i++){
       xold[i] = x[i];
     }

      promv (a1,a2,a3,a4,a5,x,r,nx,ny,ne);/* R �� AX ���� */
      for (i=1;i<=ne;i++){/* r_{0}��p_{0}(�����l)�C������ s=r_{0} �̐ݒ� */
        r[i] = b[i] - r[i]; p[i] = r[i]; s[i] = r[i];
      }

      /* �J��Ԃ��v�Z */
      for (j=1;j<=NITR;j++){
        sr1=provv (s,r,ne);/* ( s, r_{k} ) �̌v�Z => SR1 */
        promv (a1,a2,a3,a4,a5,p,ap,nx,ny,ne);/* A p_{k} �̌v�Z => AP(NE) */
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
        promv (a1,a2,a3,a4,a5,t,at,nx,ny,ne);/* A t_{k} �̌v�Z => AT(NE) */
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
