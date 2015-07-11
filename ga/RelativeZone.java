/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

/**
 *
 * @author mar
 */
//import java.util.ArrayList;
//import java.util.Random;
//import java.util.Enumeration;
public class RelativeZone{
    private FCCLattice3D fccLatt;
    public RelativeZone(){
        fccLatt=new FCCLattice3D();
    }
    public String encode(AminoAcid[] aa){

        String result="";
        int [][] baseMatrix=new int[3][3];

        baseMatrix[0][0]=1;baseMatrix[0][1]=0;baseMatrix[0][2]=0;
        baseMatrix[1][0]=0;baseMatrix[1][1]=1;baseMatrix[1][2]=0;
        baseMatrix[2][0]=0;baseMatrix[2][1]=0;baseMatrix[2][2]=1;
        // re calculate the chain to see if FCC decoding was correct
        int[] abs=new int[3];
        
        int i=0;
        for(i=0;i<Global.seqLength-1;i++){
                abs[0]=aa[i+1].x-aa[i].x;
                abs[1]=aa[i+1].y-aa[i].y;
                abs[2]=aa[i+1].z-aa[i].z;
                int idx=fccLatt.getRelativeVector(abs, baseMatrix);
                if(idx<0||idx>11){
                    Global.SomethingWrong=true;
                    return "";
                }

                char ch=Global.baseDir[idx];
                fccLatt.upadteBaseMatrixInverse(idx, baseMatrix);
                result+=ch;

        }

        return result;
    }

    public Conformation decode(String encodedString){
	int seqLen=Global.seqLength;
        AminoAcid[] aa=new AminoAcid[seqLen];
        int x=0,y=0,z=0;
        Node3d currnode=new Node3d(x, y, z);        
        char[] dir=encodedString.toCharArray();
        Global.occupiedNode.clear();

        FCCLattice3D fcclatt=new FCCLattice3D();

        int[] currPoint={0,0,0};
        aa[0]=new AminoAcid(currPoint[0],currPoint[1],currPoint[2]);
        Global.occupiedNode.add(new Node3d(currPoint[0],currPoint[1],currPoint[2]));
        int [][] baseMatrix={{1,0,0 },{0,1,0 },{0,0,1}};
        for(int i=0;i<seqLen-1;i++)
        {
            fcclatt.getNextPoint(dir[i]-97, baseMatrix,currPoint);

    //        conf.lattice->getNextPoint(conf.basisIdx[i]-'A',baseMatrix,currPoint); // get Next point into tempAxes
            fcclatt.updateBaseMatrix(dir[i]-97,baseMatrix);
            currnode=new Node3d(currPoint[0],currPoint[1],currPoint[2]);
            if (!Common.isFree(currnode,Global.occupiedNode)){
                    return null;
                }else{
                   aa[i+1]=new AminoAcid(currPoint[0],currPoint[1],currPoint[2]);
                    Global.occupiedNode.add(currnode);
                }
          }
//          return aa;
        String Seq=encode(aa);
//        int fitness=Common.calculateFitness(aa,0);
        int fitness[]={0,0};
        Common.calculateFitness(aa,fitness);
        Conformation c = new Conformation(aa, fitness, Seq);

        return c;
    }

}
