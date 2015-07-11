/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ga;

//import java.util.ArrayList;
//import java.util.ArrayList;

/**
 *
 * @author mar
 */
public class FCCLattice3D {

   
    public Global.NeighborData[] fccNeighbor;
    public char[] basisIdx;
//    public char[] baseDir=new char[12];
//    public int[][] baseVec;
//    public int[][][] baseMat;
//    public int [][][] inverMat;
    
    FCCLattice3D() {
        fccNeighbor = new Global.NeighborData[12];
        basisIdx=Global.twelveDirections.toCharArray();
        //String[] name = {"FL-0", "LU-1", "FU-2", "BL-3", "RU-4", "BU-5", "FR-6", "LD-7", "FD-8", "BR-9", "RD-10", "BD-11"};
        char[] dir = "abcdefghijkl".toCharArray();
        Global.baseDir=dir;
       int[][] vec = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {-1, 0, 1},
                {1, -1, 0}, {0, 1, -1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}};
       Global.baseVec=vec;

        int[][][] mat = {
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
            {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},
            {{1, 0, 0}, {0, 0, -1}, {0, 1, 0}},
            {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
            {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}},
            {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}},
            {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
            {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}},
            {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
            {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
            {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}},
            {{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}},};
        Global.baseMat=mat;

        int[][][] invmat = {
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
            {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}},
            {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
            {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
            {{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}},
            {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}},
            {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
            {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},
            {{1, 0, 0}, {0, 0, -1}, {0, 1, 0}},
            {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
            {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}},
            {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}}
        };
        Global.invrMat=invmat;
        for (int i = 0; i < 12; i++) {
            Node3d[] m = new Node3d[3];
            Node3d[] invm = new Node3d[3];

            int x, y, z;
            x = vec[i][0];
            y = vec[i][1];
            z = vec[i][2];
            Node3d v= new Node3d(x, y, z);
            //char c =dir[i])
            for (int j = 0; j < 3; j++) {
                x = mat[i][j][0];
                y = mat[i][j][1];
                z = mat[i][j][2];
                m[j] = new Node3d(x, y, z);
                x = invmat[i][j][0];
                y = invmat[i][j][1];
                z = invmat[i][j][2];
                invm[j] = new Node3d(x, y, z);
            }
            fccNeighbor[i]=new Global.NeighborData(dir[i], v, m, invm);
        }
    }

//    public Node3d getNextPoint(int idx, Node3d baseMatrix[]) {
//        // from basematrix and point calculate new point and save into point
//        //Node3d temp;
//        int x = 0, y = 0, z = 0;
//        for (int i = 0; i < 3; i++) {
//            x = x + fccNeighbor[idx].vec.x * baseMatrix[i].x;
//            y = y + fccNeighbor[idx].vec.y * baseMatrix[i].y;
//            z = z + fccNeighbor[idx].vec.z * baseMatrix[i].z;
//        }
//        return new Node3d(x, y, z);
//    }
public void getNextPoint(int idx, int [][] baseMatrix,int [] point) {

    // from basematrix and point calculate new point and save into point
    int[] temp={0,0,0};
    int i,j;
    for(i=0;i<3;i++)
    {
        temp[i]=0;
        for(j=0;j<3;j++)
        {
            temp[i]=temp[i]+Global.baseVec[idx][j]*baseMatrix[i][j];
        }
        point[i]=point[i]+temp[i];
    }
}

    public int getFCCbasis(int dirX, int dirY, int dirZ) {

        for (int i=0;i<12;i++){
            if(Global.baseVec[i][0]==dirX && Global.baseVec[i][1]==dirY && Global.baseVec[i][2]==dirZ ){
                return i;
            }
        }
        return -1;
    }

    public int getRelativeVector(Node3d abs, Node3d[] BaseMatrix) {
        int x = 0, y = 0, z = 0;
        for (int i = 0; i < 3; i++) {
            x = x + BaseMatrix[i].x * abs.x;
            y = y + BaseMatrix[i].y * abs.y;
            z = z + BaseMatrix[i].z * abs.z;
        }
        return getFCCbasis(x, y, z);
    }

     public int getRelativeVector(int[] abs, int[][] BaseMatrix) {
        int[] temp=new int[3];
        int i,j;
        for(i=0;i<3;i++)
        {
            temp[i]=0;
            for(j=0;j<3;j++)
            {
                temp[i]=temp[i]+BaseMatrix[i][j]*abs[j];
            }
        }
        return getFCCbasis(temp[0],temp[1],temp[2]);
    }

    public Node3d[] upadteBaseMatrixInverse(int idx, Node3d[] BaseMatrix) {
        //int idx=getFCCbasis(abs[0],abs[1],abs[2]);
        //Node3d[] temp=new Node3d[3];
        //int i,j,k;
        int[] cord=new int[3];

        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                cord[j]=fccNeighbor[idx].invmat[i].x * BaseMatrix[j].x
                        +fccNeighbor[idx].invmat[i].y * BaseMatrix[j].y
                        +fccNeighbor[idx].invmat[i].z * BaseMatrix[j].z;
            }
            BaseMatrix[i] = new Node3d(cord[0],cord[1],cord[2]);
        }
        return BaseMatrix;
    }
    public void upadteBaseMatrixInverse(int idx, int[][] BaseMatrix) {
        //int idx=getFCCbasis(abs[0],abs[1],abs[2]);
    int[][] temp=new int[3][3];
    int i,j,k;
    for(i=0;i<3;i++)
    {

        for(j=0;j<3;j++)
        {
            temp[i][j]=0;
            for(k=0;k<3;k++)
            {
                 temp[i][j]=temp[i][j]+Global.invrMat[idx][i][k]*BaseMatrix[k][j];
            }
        }
    }
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            BaseMatrix[i][j]=temp[i][j];
        }
    }
    }
    public void updateBaseMatrix(int idx, Node3d[] BaseMatrix) {

        for (int i = 0; i < 3; i++) {
            int x = 0, y = 0, z = 0;
            x = x + fccNeighbor[idx].mat[i].x * BaseMatrix[i].x;
            y = y + fccNeighbor[idx].mat[i].y * BaseMatrix[i].y;
            z = y + fccNeighbor[idx].mat[i].z * BaseMatrix[i].z;

            BaseMatrix[i] = new Node3d(x, y, z);
        }
    }
    public void updateBaseMatrix(int idx, int[][] BaseMatrix) {
// from basaematrix and the move matrix calculate new base matrix
        int [][] temp=new int[3][3];
        int i,j,k;
        for(i=0;i<3;i++)
        {

            for(j=0;j<3;j++)
            {
                temp[i][j]=0;
                for(k=0;k<3;k++)
                {
                    temp[i][j]=temp[i][j]+BaseMatrix[i][k]*Global.baseMat[idx][k][j];
                }
            }
        }
        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                BaseMatrix[i][j]=temp[i][j];
            }
        }
    }
//    public int getAbsoluteVectorByDim(int idx,int dim)
//    {
//            return fccNeighborData[idx].vec[dim];
//    }

    public void getAbsoluteDirection(int idx,int[]abs)
    {
        //return fccNeighbor[idx].vec;
        for(int i=0;i<3;i++){
            abs[i]=Global.baseVec[idx][i];
        }
//        abs[1]=Global.baseVec[idx][0];
//        abs[2]=Global.baseVec[idx][0];
    }

//    public static void main(String[] args) {
//       FCCLattice3D f=new FCCLattice3D();
//       for (int i=0;i<12;i++){
//            System.out.println(Global.baseDir[i]);
//        }
//       for (int i=0;i<12;i++){
//           String s="";
//           for (int j=0;j<3;j++){
//            s=s+Global.baseVec[i][j];
//           }
//           System.out.println(s);
//
//        }
//       System.out.println("===================");
//        for (int i=0;i<12;i++){
//
//           for (int j=0;j<3;j++){
//               String s="";
//               for (int k=0;k<3;k++){
//                    s=s+Global.baseMat[i][j][k];
//               }
//               System.out.println(s);
//           }
//           System.out.println("===================");
//        }
//       for (int i=0;i<12;i++){
//
//           for (int j=0;j<3;j++){
//               String s="";
//               for (int k=0;k<3;k++){
//                    s=s+Global.invrMat[i][j][k];
//               }
//               System.out.println(s);
//           }
//           System.out.println("===================");
//        }
//    }



}
