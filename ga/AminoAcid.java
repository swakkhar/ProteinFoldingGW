package ga;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Administrator
 */
public class AminoAcid {
    public int x;
    public int y;
    public int z;
    public int[] operation=new int[5];
    AminoAcid()
    {
        int[] op={-1,-1,-1,-1,-1};
        this.x=0;
        this.y=0;
        this.z=0;
        this.operation=op;
    }
    AminoAcid(int XCoordinate,int YCoordinate,int ZCoordinate)
    {
        int[] op={-1,-1,-1,-1,-1};
        this.x=XCoordinate;
        this.y=YCoordinate;
        this.z=ZCoordinate;
        int _i=0;
        for (int el:op)
            this.operation[_i++]=el;
    }
    AminoAcid(int XCoordinate,int YCoordinate,int ZCoordinate,int[] _operation)
    {
        this.operation=new int[_operation.length];
        this.x=XCoordinate;
        this.y=YCoordinate;
        this.z=ZCoordinate;
        int _i=0;
        for (int el:_operation)
            this.operation[_i++]=el;
    }
    AminoAcid(AminoAcid aa)
    {
        this.operation=new int[aa.operation.length];
        this.x=aa.x;
        this.y=aa.y;
        this.z=aa.z;
        int _i=0;
        for (int el:aa.operation)
            this.operation[_i++]=el;
    }

    public void copyOperation(int op[]){
        int i=0;
        for (int j:this.operation){
            op[i++]=j;
        }
    }


//    	public int x;
//        public int y;
//        public int z;
//        public char type;
//
//    	AminoAcid()
//        {
//            this.x=0;
//            this.y=0;
//            this.z=0;
//            this.type='0';
//        }
//         AminoAcid(int XCoordinate,int YCoordinate,int ZCoordinate,char Type)
//        {
//            this.x=XCoordinate;
//            this.y=YCoordinate;
//            this.z=ZCoordinate;
//            this.type=Type;
//        }
//        AminoAcid(AminoAcid aa)
//        {
//            this.x=aa.x;
//            this.y=aa.y;
//            this.z=aa.z;
//            this.type=aa.type;
//        }
}