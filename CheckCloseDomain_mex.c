/**************************************
*
*
*       Check if the domain is closed
*
*
*****************************************/

#include <math.h>
#include <mex.h>

static void CheckCloseDomain_mex(double x11, double x12, double x13,
                                double x21, double x22, double x23,
                                double x31, double x32, double x33,
                                double x41, double x42, double x43,
                                double SFT_plane1_1, double SFT_plane1_2, double SFT_plane1_3, double SFT_plane1_4, double SFT_plane1_5, double SFT_plane1_6, double SFT_plane1_7, double SFT_plane1_8, double SFT_plane1_9, double SFT_plane1_10, double SFT_plane1_11, double SFT_plane1_12,
                                double SFT_plane2_1, double SFT_plane2_2, double SFT_plane2_3, double SFT_plane2_4, double SFT_plane2_5, double SFT_plane2_6, double SFT_plane2_7, double SFT_plane2_8, double SFT_plane2_9, double SFT_plane2_10, double SFT_plane2_11, double SFT_plane2_12,
                                double SFT_plane3_1, double SFT_plane3_2, double SFT_plane3_3, double SFT_plane3_4, double SFT_plane3_5, double SFT_plane3_6, double SFT_plane3_7, double SFT_plane3_8, double SFT_plane3_9, double SFT_plane3_10, double SFT_plane3_11, double SFT_plane3_12,
                                double SFT_plane4_1, double SFT_plane4_2, double SFT_plane4_3, double SFT_plane4_4, double SFT_plane4_5, double SFT_plane4_6, double SFT_plane4_7, double SFT_plane4_8, double SFT_plane4_9, double SFT_plane4_10, double SFT_plane4_11, double SFT_plane4_12,
                                int *check);
                                
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])

{
    double *x1, *x2, *x3, *x4;
    double *SFT_plane1, *SFT_plane2, *SFT_plane3, *SFT_plane4;
    double x11, x12, x13, x21, x22, x23, x31, x32, x33, x41, x42, x43;
    int doSFT;
    double SFT_plane1_1, SFT_plane1_2, SFT_plane1_3, SFT_plane1_4, SFT_plane1_5, SFT_plane1_6, SFT_plane1_7, SFT_plane1_8, SFT_plane1_9, SFT_plane1_10, SFT_plane1_11, SFT_plane1_12;
    double SFT_plane2_1, SFT_plane2_2, SFT_plane2_3, SFT_plane2_4, SFT_plane2_5, SFT_plane2_6, SFT_plane2_7, SFT_plane2_8, SFT_plane2_9, SFT_plane2_10, SFT_plane2_11, SFT_plane2_12;
    double SFT_plane3_1, SFT_plane3_2, SFT_plane3_3, SFT_plane3_4, SFT_plane3_5, SFT_plane3_6, SFT_plane3_7, SFT_plane3_8, SFT_plane3_9, SFT_plane3_10, SFT_plane3_11, SFT_plane3_12;
    double SFT_plane4_1, SFT_plane4_2, SFT_plane4_3, SFT_plane4_4, SFT_plane4_5, SFT_plane4_6, SFT_plane4_7, SFT_plane4_8, SFT_plane4_9, SFT_plane4_10, SFT_plane4_11, SFT_plane4_12;
    int check=0;
    
    /* Check for proper number of input arguments. */
    if (nrhs!=8) mexErrMsgTxt("8 inputs required!");
    if (nlhs!=1) mexErrMsgTxt("1 outputs required!");
 
    
    x1=mxGetPr(prhs[0]);
    x2=mxGetPr(prhs[1]);
    x3=mxGetPr(prhs[2]);
    x4=mxGetPr(prhs[3]);
    
    x11=*(x1);
    x12=*(x1+1);
    x13=*(x1+2);
    x21=*(x2);
    x22=*(x2+1);
    x23=*(x2+2);
    x31=*(x3);
    x32=*(x3+1);
    x33=*(x3+2);
    x41=*(x4);
    x42=*(x4+1);
    x43=*(x4+2);
    
    SFT_plane1=mxGetPr(prhs[4]);
    SFT_plane2=mxGetPr(prhs[5]);
    SFT_plane3=mxGetPr(prhs[6]);
    SFT_plane4=mxGetPr(prhs[7]);
    
    SFT_plane1_1=*(SFT_plane1);
    SFT_plane1_2=*(SFT_plane1+1);
    SFT_plane1_3=*(SFT_plane1+2);
    SFT_plane1_4=*(SFT_plane1+3);
    SFT_plane1_5=*(SFT_plane1+4);
    SFT_plane1_6=*(SFT_plane1+5);
    SFT_plane1_7=*(SFT_plane1+6);
    SFT_plane1_8=*(SFT_plane1+7);
    SFT_plane1_9=*(SFT_plane1+8);
    SFT_plane1_10=*(SFT_plane1+9);
    SFT_plane1_11=*(SFT_plane1+10);
    SFT_plane1_12=*(SFT_plane1+11);
    
    SFT_plane2_1=*(SFT_plane2);
    SFT_plane2_2=*(SFT_plane2+1);
    SFT_plane2_3=*(SFT_plane2+2);
    SFT_plane2_4=*(SFT_plane2+3);
    SFT_plane2_5=*(SFT_plane2+4);
    SFT_plane2_6=*(SFT_plane2+5);
    SFT_plane2_7=*(SFT_plane2+6);
    SFT_plane2_8=*(SFT_plane2+7);
    SFT_plane2_9=*(SFT_plane2+8);
    SFT_plane2_10=*(SFT_plane2+9);
    SFT_plane2_11=*(SFT_plane2+10);
    SFT_plane2_12=*(SFT_plane2+11);
    
    SFT_plane3_1=*(SFT_plane3);
    SFT_plane3_2=*(SFT_plane3+1);
    SFT_plane3_3=*(SFT_plane3+2);
    SFT_plane3_4=*(SFT_plane3+3);
    SFT_plane3_5=*(SFT_plane3+4);
    SFT_plane3_6=*(SFT_plane3+5);
    SFT_plane3_7=*(SFT_plane3+6);
    SFT_plane3_8=*(SFT_plane3+7);
    SFT_plane3_9=*(SFT_plane3+8);
    SFT_plane3_10=*(SFT_plane3+9);
    SFT_plane3_11=*(SFT_plane3+10);
    SFT_plane3_12=*(SFT_plane3+11);
    
    SFT_plane4_1=*(SFT_plane4);
    SFT_plane4_2=*(SFT_plane4+1);
    SFT_plane4_3=*(SFT_plane4+2);
    SFT_plane4_4=*(SFT_plane4+3);
    SFT_plane4_5=*(SFT_plane4+4);
    SFT_plane4_6=*(SFT_plane4+5);
    SFT_plane4_7=*(SFT_plane4+6);
    SFT_plane4_8=*(SFT_plane4+7);
    SFT_plane4_9=*(SFT_plane4+8);
    SFT_plane4_10=*(SFT_plane4+9);
    SFT_plane4_11=*(SFT_plane4+10);
    SFT_plane4_12=*(SFT_plane4+11);
    
    CheckCloseDomain_mex(x11,x12,x13,
    x21,x22,x23,
    x31,x32,x33,
    x41,x42,x43,
    SFT_plane1_1,SFT_plane1_2,SFT_plane1_3,SFT_plane1_4,SFT_plane1_5,SFT_plane1_6,SFT_plane1_7,SFT_plane1_8,SFT_plane1_9,SFT_plane1_10,SFT_plane1_11,SFT_plane1_12,
    SFT_plane2_1,SFT_plane2_2,SFT_plane2_3,SFT_plane2_4,SFT_plane2_5,SFT_plane2_6,SFT_plane2_7,SFT_plane2_8,SFT_plane2_9,SFT_plane2_10,SFT_plane2_11,SFT_plane2_12,
    SFT_plane3_1,SFT_plane3_2,SFT_plane3_3,SFT_plane3_4,SFT_plane3_5,SFT_plane3_6,SFT_plane3_7,SFT_plane3_8,SFT_plane3_9,SFT_plane3_10,SFT_plane3_11,SFT_plane3_12,
    SFT_plane4_1,SFT_plane4_2,SFT_plane4_3,SFT_plane4_4,SFT_plane4_5,SFT_plane4_6,SFT_plane4_7,SFT_plane4_8,SFT_plane4_9,SFT_plane4_10,SFT_plane4_11,SFT_plane4_12,
    &check);
    
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=check;
        
}


        
        
        /**
         *
         * Check if the Domain is closed
         *
         **/
        
static void CheckCloseDomain_mex(double x11, double x12, double x13,
                                double x21, double x22, double x23,
                                double x31, double x32, double x33,
                                double x41, double x42, double x43,
                                double SFT_plane1_1, double SFT_plane1_2, double SFT_plane1_3, double SFT_plane1_4, double SFT_plane1_5, double SFT_plane1_6, double SFT_plane1_7, double SFT_plane1_8, double SFT_plane1_9, double SFT_plane1_10, double SFT_plane1_11, double SFT_plane1_12,
                                double SFT_plane2_1, double SFT_plane2_2, double SFT_plane2_3, double SFT_plane2_4, double SFT_plane2_5, double SFT_plane2_6, double SFT_plane2_7, double SFT_plane2_8, double SFT_plane2_9, double SFT_plane2_10, double SFT_plane2_11, double SFT_plane2_12,
                                double SFT_plane3_1, double SFT_plane3_2, double SFT_plane3_3, double SFT_plane3_4, double SFT_plane3_5, double SFT_plane3_6, double SFT_plane3_7, double SFT_plane3_8, double SFT_plane3_9, double SFT_plane3_10, double SFT_plane3_11, double SFT_plane3_12,
                                double SFT_plane4_1, double SFT_plane4_2, double SFT_plane4_3, double SFT_plane4_4, double SFT_plane4_5, double SFT_plane4_6, double SFT_plane4_7, double SFT_plane4_8, double SFT_plane4_9, double SFT_plane4_10, double SFT_plane4_11, double SFT_plane4_12,
                                int *check)
                                
{
    double seg[4][6];
    double SFT_plane[4][12];
    double seg_unit_vector[3];
    double plane_normal[3];
    double norm_plane=0;
    double norm_vec=0;
    double cut_point[3];
    double norm_cut_point=0;
    double tol=0, delta=0;
    double a[3], b[3], c[3];
    double ab_vec[3], ac_vec[3], bc_vec[3];
    double ab_vec_norm=0, ac_vec_norm=0, bc_vec_norm=0;
    double vec1[3], vec2[3], vec3[3];
    double x1[3], x2[3];
    double min_x=0, max_x=0;
    double norm_x1=0, norm_x2=0;
    double norm_vec1=0, norm_vec2=0, norm_vec3=0;
    double theta1=0, theta2=0, theta3=0;
    double dotvec1vec2=0, dotvec2vec3=0, dotvec3vec1=0;
    int j=0, k=0;
    int point_inside=0;
    
    tol=M_PI/100;
    delta=0.1;
        
    seg[0][0]=x11;
    seg[0][1]=x12;
    seg[0][2]=x13;
    seg[0][3]=x31;
    seg[0][4]=x32;
    seg[0][5]=x33;
    
    seg[1][0]=x11;
    seg[1][1]=x12;
    seg[1][2]=x13;
    seg[1][3]=x41;
    seg[1][4]=x42;
    seg[1][5]=x43;
    
    seg[2][0]=x21;
    seg[2][1]=x22;
    seg[2][2]=x23;
    seg[2][3]=x31;
    seg[2][4]=x32;
    seg[2][5]=x33;
    
    seg[3][0]=x21;
    seg[3][1]=x22;
    seg[3][2]=x23;
    seg[3][3]=x41;
    seg[3][4]=x42;
    seg[3][5]=x43;
    
    SFT_plane[0][0]=SFT_plane1_1;
    SFT_plane[0][1]=SFT_plane1_2;
    SFT_plane[0][2]=SFT_plane1_3;
    SFT_plane[0][3]=SFT_plane1_4;
    SFT_plane[0][4]=SFT_plane1_5;
    SFT_plane[0][5]=SFT_plane1_6;
    SFT_plane[0][6]=SFT_plane1_7;
    SFT_plane[0][7]=SFT_plane1_8;
    SFT_plane[0][8]=SFT_plane1_9;
    SFT_plane[0][9]=SFT_plane1_10;
    SFT_plane[0][10]=SFT_plane1_11;
    SFT_plane[0][11]=SFT_plane1_12;
    
    SFT_plane[1][0]=SFT_plane2_1;
    SFT_plane[1][1]=SFT_plane2_2;
    SFT_plane[1][2]=SFT_plane2_3;
    SFT_plane[1][3]=SFT_plane2_4;
    SFT_plane[1][4]=SFT_plane2_5;
    SFT_plane[1][5]=SFT_plane2_6;
    SFT_plane[1][6]=SFT_plane2_7;
    SFT_plane[1][7]=SFT_plane2_8;
    SFT_plane[1][8]=SFT_plane2_9;
    SFT_plane[1][9]=SFT_plane2_10;
    SFT_plane[1][10]=SFT_plane2_11;
    SFT_plane[1][11]=SFT_plane2_12;
        
    SFT_plane[2][0]=SFT_plane3_1;
    SFT_plane[2][1]=SFT_plane3_2;
    SFT_plane[2][2]=SFT_plane3_3;
    SFT_plane[2][3]=SFT_plane3_4;
    SFT_plane[2][4]=SFT_plane3_5;
    SFT_plane[2][5]=SFT_plane3_6;
    SFT_plane[2][6]=SFT_plane3_7;
    SFT_plane[2][7]=SFT_plane3_8;
    SFT_plane[2][8]=SFT_plane3_9;
    SFT_plane[2][9]=SFT_plane3_10;
    SFT_plane[2][10]=SFT_plane3_11;
    SFT_plane[2][11]=SFT_plane3_12;    
        
    SFT_plane[3][0]=SFT_plane4_1;
    SFT_plane[3][1]=SFT_plane4_2;
    SFT_plane[3][2]=SFT_plane4_3;
    SFT_plane[3][3]=SFT_plane4_4;
    SFT_plane[3][4]=SFT_plane4_5;
    SFT_plane[3][5]=SFT_plane4_6;
    SFT_plane[3][6]=SFT_plane4_7;
    SFT_plane[3][7]=SFT_plane4_8;
    SFT_plane[3][8]=SFT_plane4_9;
    SFT_plane[3][9]=SFT_plane4_10;
    SFT_plane[3][10]=SFT_plane4_11;
    SFT_plane[3][11]=SFT_plane4_12;    
        
    for(j=0;j<4;j++){  /* Each plane of the SFT */
        
        norm_plane=sqrt(pow(SFT_plane[j][0],2)+pow(SFT_plane[j][1],2)+pow(SFT_plane[j][2],2));
        plane_normal[0]=SFT_plane[j][0]/norm_plane;
        plane_normal[1]=SFT_plane[j][1]/norm_plane;
        plane_normal[2]=SFT_plane[j][2]/norm_plane;
        
        a[0]=SFT_plane[j][3];
        a[1]=SFT_plane[j][4];
        a[2]=SFT_plane[j][5];
        
        b[0]=SFT_plane[j][6];
        b[1]=SFT_plane[j][7];
        b[2]=SFT_plane[j][8];
        
        c[0]=SFT_plane[j][9];
        c[1]=SFT_plane[j][10];
        c[2]=SFT_plane[j][11];
        
        ab_vec_norm=sqrt(pow((b[0]-a[0]),2) + pow((b[1]-a[1]),2) + pow((b[2]-a[2]),2));
        ac_vec_norm=sqrt(pow((c[0]-a[0]),2) + pow((c[1]-a[1]),2) + pow((c[2]-a[2]),2));
        bc_vec_norm=sqrt(pow((c[0]-b[0]),2) + pow((c[1]-b[1]),2) + pow((c[2]-b[2]),2));
        
        ab_vec[0]=(b[0]-a[0])/ab_vec_norm;
        ab_vec[1]=(b[1]-a[1])/ab_vec_norm;
        ab_vec[2]=(b[2]-a[2])/ab_vec_norm;
        
        ac_vec[0]=(c[0]-a[0])/ac_vec_norm;
        ac_vec[1]=(c[1]-a[1])/ac_vec_norm;
        ac_vec[2]=(c[2]-a[2])/ac_vec_norm;
        
        bc_vec[0]=(c[0]-b[0])/bc_vec_norm;
        bc_vec[1]=(c[1]-b[1])/bc_vec_norm;
        bc_vec[2]=(c[2]-b[2])/bc_vec_norm;
        
        a[0]=SFT_plane[j][3] + 2*delta*(ab_vec[0]+ac_vec[0]) - delta*plane_normal[0];
        a[1]=SFT_plane[j][4] + 2*delta*(ab_vec[1]+ac_vec[1]) - delta*plane_normal[1];
        a[2]=SFT_plane[j][5] + 2*delta*(ab_vec[2]+ac_vec[2]) - delta*plane_normal[2];
        
        b[0]=SFT_plane[j][6] + 2*delta*(bc_vec[0]-ab_vec[0]) - delta*plane_normal[0];
        b[1]=SFT_plane[j][7] + 2*delta*(bc_vec[1]-ab_vec[1]) - delta*plane_normal[1];
        b[2]=SFT_plane[j][8] + 2*delta*(bc_vec[2]-ab_vec[2]) - delta*plane_normal[2];
        
        c[0]=SFT_plane[j][9] + 2*delta*(-bc_vec[0]-ac_vec[0]) - delta*plane_normal[0];
        c[1]=SFT_plane[j][10] + 2*delta*(-bc_vec[1]-ac_vec[1]) - delta*plane_normal[1];
        c[2]=SFT_plane[j][11] + 2*delta*(-bc_vec[2]-ac_vec[2]) - delta*plane_normal[2];
        
        for(k=0;k<4;k++){  /* Each segment formed by each pair of nodes for each segment */
            
            if((seg[k][0]!=seg[k][3])||(seg[k][1]!=seg[k][4])||(seg[k][2]!=seg[k][5])){
                
                point_inside=0;
                x1[0]=seg[k][0];
                x1[1]=seg[k][1];
                x1[2]=seg[k][2];
                norm_x1=sqrt(pow(x1[0],2) + pow(x1[1],2) + pow(x1[2],2));
                
                x2[0]=seg[k][3];
                x2[1]=seg[k][4];
                x2[2]=seg[k][5];
                norm_x2=sqrt(pow(x2[0],2) + pow(x2[1],2) + pow(x2[2],2));
                
                norm_vec=sqrt(pow((seg[k][0]-seg[k][3]),2)+pow((seg[k][1]-seg[k][4]),2)+pow((seg[k][2]-seg[k][5]),2));
                seg_unit_vector[0]=(seg[k][0]-seg[k][3])/norm_vec;
                seg_unit_vector[1]=(seg[k][1]-seg[k][4])/norm_vec;
                seg_unit_vector[2]=(seg[k][2]-seg[k][5])/norm_vec;
                
                cut_point[2]=((plane_normal[0]*(seg_unit_vector[0]/seg_unit_vector[2]) + plane_normal[1]*(seg_unit_vector[1]/seg_unit_vector[2]))*seg[k][2] + plane_normal[2]*a[2] - plane_normal[0]*(seg[k][0] - a[0]) - plane_normal[1]*(seg[k][1] - a[1]))/(plane_normal[0]*(seg_unit_vector[0]/seg_unit_vector[2]) + plane_normal[1]*(seg_unit_vector[1]/seg_unit_vector[2]) + plane_normal[2]);
                cut_point[0]=(seg_unit_vector[0]/seg_unit_vector[2])*(cut_point[2] - seg[k][2]) + seg[k][0];
                cut_point[1]=(seg_unit_vector[1]/seg_unit_vector[2])*(cut_point[2] - seg[k][2]) + seg[k][1];
                norm_cut_point=sqrt(pow(cut_point[0],2) + pow(cut_point[1],2) + pow(cut_point[2],2));
                
                norm_vec1=sqrt(pow((a[0] - cut_point[0]),2)+pow((a[1] - cut_point[1]),2) + pow((a[2] - cut_point[2]),2));
                vec1[0]=(a[0] - cut_point[0])/norm_vec1;
                vec1[1]=(a[1] - cut_point[1])/norm_vec1;
                vec1[2]=(a[2] - cut_point[2])/norm_vec1;
                
                norm_vec2=sqrt(pow((b[0] - cut_point[0]),2)+pow((b[1] - cut_point[1]),2) + pow((b[2] - cut_point[2]),2));
                vec2[0]=(b[0] - cut_point[0])/norm_vec2;
                vec2[1]=(b[1] - cut_point[1])/norm_vec2;
                vec2[2]=(b[2] - cut_point[2])/norm_vec2;
                
                norm_vec3=sqrt(pow((c[0] - cut_point[0]),2)+pow((c[1] - cut_point[1]),2) + pow((c[2] - cut_point[2]),2));
                vec3[0]=(c[0] - cut_point[0])/norm_vec3;
                vec3[1]=(c[1] - cut_point[1])/norm_vec3;
                vec3[2]=(c[2] - cut_point[2])/norm_vec3;
                
                dotvec1vec2=vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
                dotvec2vec3=vec2[0]*vec3[0] + vec2[1]*vec3[1] + vec2[2]*vec3[2];
                dotvec3vec1=vec3[0]*vec1[0] + vec3[1]*vec1[1] + vec3[2]*vec1[2];
                
                theta1=acos(dotvec1vec2);
                theta2=acos(dotvec2vec3);
                theta3=acos(dotvec3vec1);
                
                if(((2*M_PI)-(theta1+theta2+theta3))<tol){
                    point_inside=1;
                }
                else{
                    point_inside=0;
                    continue;
                }
                
                min_x=norm_x1;
                max_x=norm_x2;
                if(min_x>max_x){
                    min_x=norm_x2;
                    max_x=norm_x1;
                }
                
                if((point_inside)&&(norm_cut_point>min_x)&&(norm_cut_point<max_x)){
                    *check=1;
                    return;
                }
                
            }
            
        }
        
    }
    
}
