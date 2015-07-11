/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

import java.util.ArrayList;
//import java.util.Collections;
//import java.util.Random;
//import javax.swing.JFrame;
//
//import java.awt.Frame;
//import java.awt.event.*;
//import java.awt.GraphicsConfiguration;
//import com.sun.j3d.utils.applet.MainFrame;
//import java.awt.Color;
//import sun.misc.ASCIICaseInsensitiveComparator;
public class FCCDescriptor {

private static ArrayList<Node3d> relativePoints=new ArrayList<Node3d>();

private static ArrayList<Node3d> A=new ArrayList<Node3d>();
private static ArrayList<Node3d> B=new ArrayList<Node3d>();
private static ArrayList<Node3d> C=new ArrayList<Node3d>();
private static ArrayList<Node3d> D=new ArrayList<Node3d>();
private static ArrayList<Node3d> E=new ArrayList<Node3d>();
private static ArrayList<Node3d> F=new ArrayList<Node3d>();
private static ArrayList<Node3d> G=new ArrayList<Node3d>();
private static ArrayList<Node3d> H=new ArrayList<Node3d>();
private static ArrayList<Node3d> I=new ArrayList<Node3d>();
private static ArrayList<Node3d> J=new ArrayList<Node3d>();
private static ArrayList<Node3d> K=new ArrayList<Node3d>();
private static ArrayList<Node3d> L=new ArrayList<Node3d>();

private static int sine(int angle)
{
	if(angle==90)
		return 1;
	else if(angle==0)
		return 0;
	else if(angle==-90)
		return -1;
	else if(angle==-180)
		return 0;
	else if(angle==180)
		return 0;
        else
            return -11;
}
private static int cosine(int angle)
{
	if(angle==90)
		return 0;
	else if(angle==0)
		return 1;
	else if(angle==-90)
		return 0;
	else if(angle==-180)
		return -1;
	else if(angle==180)
		return -1;
        else
            return -11;
}
private static Node3d xrotate(Node3d p1, int theta)
{
    Node3d p2;
    int x,y,z;
    int tempx,tempy,tempz;
    tempx=p1.x;tempy=p1.y;tempz=p1.z;

    x=tempx;
    y=tempy*cosine(theta)-tempz*sine(theta);
    z=tempy*sine(theta)+tempz*cosine(theta);
    p2=new Node3d(x,y,z);
    return p2;
}
private static Node3d yrotate(Node3d p1, int theta)
{
    Node3d p2;
    int x,y,z;
    int tempx,tempy,tempz;
    tempx=p1.x;tempy=p1.y;tempz=p1.z;

   
    y=tempy;
    z=tempz*cosine(theta)-tempx*sine(theta);
    x=tempz*sine(theta)+tempx*cosine(theta);
    p2=new Node3d(x,y,z);
    return p2;
}

private static Node3d zrotate(Node3d p1, int theta)
{
    Node3d p2;
    int x,y,z;
    int tempx,tempy,tempz;
    tempx=p1.x;tempy=p1.y;tempz=p1.z;

    z=tempz;
    x=tempx*cosine(theta)-tempy*sine(theta);
    y=tempx*sine(theta)+tempy*cosine(theta);
    
    p2=new Node3d(x,y,z);
    return p2;
}
private static void Print(Node3d point)
{
    relativePoints.add(point);
    //System.out.println("{"+point.x+","+point.y+","+point.z+"}");
}
private static void generatePoints() {
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
    relativePoints.clear();
    Node3d[] p=new Node3d[12];
   //Node3d p=new Node3d(0, 0, 0);

    p[0]=new Node3d(+1,+0,-1);
    p[1]=new Node3d(-1,+0,+1);
    p[2]=new Node3d(+1,+0,+1);
    p[3]=new Node3d(-1,+0,-1);

    p[4]=new Node3d(+0,+1,+1);
    p[5]=new Node3d(+0,-1,-1);
    p[6]=new Node3d(+0,+1,-1);
    p[7]=new Node3d(+0,-1,+1);

    p[8]=new Node3d(+1,+1,+0);
    p[9]=new Node3d(-1,-1,+0);
    p[10]=new Node3d(-1,+1,+0);
    p[11]=new Node3d(+1,-1,+0);


//	//1
	for(int i=0;i<12;i++){
            relativePoints.add(new Node3d(p[i].x,p[i].y,p[i].z));
            //Print(p[i]);
    }

//	// 2
    for(int i=0;i<12;i++)
    {
        Node3d p1=yrotate(p[i],180);
        relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
        //Print(p[i]);
        yrotate(p[i],-180);

    }
		

    //3
    for(int i=0;i<12;i++)
    {
        Node3d p1=yrotate(p[i],90);
        relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
        //Print(p[i]);
        yrotate(p[i],-90);

    }

//4
    for(int i=0;i<12;i++)
    {
        Node3d p1=yrotate(p[i],-90);
        relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
        //Print(p[i]);
        yrotate(p[i],90);

    }


    //5
    for(int i=0;i<12;i++)
    {
        Node3d p1=zrotate(p[i],-90);
        relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
        //Print(p[i]);
        zrotate(p[i],90);
    }


// 6

for(int i=0;i<12;i++)
{
    Node3d p1=zrotate(p[i],-90);
    p1=xrotate(p1,180);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
    
//    xrotate(p[i],-180);
//    zrotate(p[i],90);
}


//7

for(int i=0;i<12;i++)
{
    Node3d p1=zrotate(p[i],90);
    p1=xrotate(p1,180);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    xrotate(p[i],-180);
//    zrotate(p[i],-90);
}
//8
for(int i=0;i<12;i++)
{
    Node3d p1=zrotate(p[i],-90);
    p1=xrotate(p1,90);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    xrotate(p[i],-90);
//    zrotate(p[i],90);
}
//9
for(int i=0;i<12;i++)
{
    Node3d p1=xrotate(p[i],-90);
    p1=zrotate(p1,-90);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    zrotate(p[i],90);
//    xrotate(p[i],90);
}


//10
for(int i=0;i<12;i++)
{
    Node3d p1=xrotate(p[i],-90);
    p1=zrotate(p1,90);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    zrotate(p[i],-90);
//    xrotate(p[i],90);
}

//11
for(int i=0;i<12;i++)
{
    Node3d p1=xrotate(p[i],-90);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    xrotate(p[i],90);
}

//12
for(int i=0;i<12;i++)
{
    Node3d p1=xrotate(p[i],-90);
    p1=zrotate(p1,180);
    relativePoints.add(new Node3d(p1.x,p1.y,p1.z));

//    zrotate(p[i],-180);
//    xrotate(p[i],90);
}

			
			
			/*for(int i=0;i<12;i++)
					Print(p[i]);
				*/
}
public static ArrayList get12x12BasisVectors(){
    generate12guide();
    ArrayList<ArrayList> myArrays=new ArrayList<ArrayList>();

    myArrays.add(A);
    myArrays.add(B);
    myArrays.add(C);
    myArrays.add(D);
    myArrays.add(E);
    myArrays.add(F);
    myArrays.add(G);
    myArrays.add(H);
    myArrays.add(I);
    myArrays.add(J);
    myArrays.add(K);
    myArrays.add(L);

//ArrayList<Node3d> A1=new ArrayList<Node3d>();
//ArrayList<Node3d> B1=new ArrayList<Node3d>();
//A1=myArrays.get(0);
//B1=myArrays.get(1);

    return myArrays;
}

private static void generate12guide() {
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
    relativePoints.clear();
    A.clear();B.clear();C.clear();D.clear();E.clear();F.clear();G.clear();H.clear();I.clear();J.clear();K.clear();L.clear();
    
    Node3d[] p=new Node3d[12];
   //Node3d p=new Node3d(0, 0, 0);

    p[0]=new Node3d(+1,+0,-1);
    p[1]=new Node3d(-1,+0,+1);
    p[2]=new Node3d(+1,+0,+1);
    p[3]=new Node3d(-1,+0,-1);

    p[4]=new Node3d(+0,+1,+1);
    p[5]=new Node3d(+0,-1,-1);
    p[6]=new Node3d(+0,+1,-1);
    p[7]=new Node3d(+0,-1,+1);

    p[8]=new Node3d(+1,+1,+0);
    p[9]=new Node3d(-1,-1,+0);
    p[10]=new Node3d(-1,+1,+0);
    p[11]=new Node3d(+1,-1,+0);

    for(int i=0;i<12;i++){
        Node3d p1=new Node3d(p[i].x,p[i].y,p[i].z);
        putAt(p1, i);

        p1=yrotate(p[i],180);
        putAt(p1, i);

        p1=yrotate(p[i],90);
        //relativePoints.add(new Node3d(p1.x,p1.y,p1.z));
        putAt(p1, i);

        p1=yrotate(p[i],-90);
        putAt(p1, i);

        p1=zrotate(p[i],-90);
        putAt(p1, i);

        p1=zrotate(p[i],-90);
        p1=xrotate(p1,180);
        putAt(p1, i);

        p1=zrotate(p[i],90);
        p1=xrotate(p1,180);
        putAt(p1, i);

        p1=zrotate(p[i],-90);
        p1=xrotate(p1,90);
        putAt(p1, i);

        p1=xrotate(p[i],-90);
        p1=zrotate(p1,-90);
        putAt(p1, i);

        p1=xrotate(p[i],-90);
        p1=zrotate(p1,90);
        putAt(p1, i);

        p1=xrotate(p[i],-90);
        putAt(p1, i);

        p1=xrotate(p[i],-90);
        p1=zrotate(p1,180);
        putAt(p1, i);
    }
}
private static void putAt(Node3d node, int idx){
    //Node3d node=new Node3d(point.x,point.y,point.z);

    switch(idx){
        case 0:
            A.add(node);
            break;
        case 1:
            B.add(node);
            break;
        case 2:
            C.add(node);
            break;
        case 3:
            D.add(node);
            break;
        case 4:
            E.add(node);
            break;
        case 5:
            F.add(node);
            break;
        case 6:
            G.add(node);
            break;
        case 7:
            H.add(node);
            break;
        case 8:
            I.add(node);
            break;
        case 9:
            J.add(node);
            break;
        case 10:
            K.add(node);
            break;
        case 11:
            L.add(node);
            break;
    }
}
private static void printAarry(ArrayList<Node3d> node){
    String s4="",s5="",s6="";
    String comma;
    for(int j=0;j<12;j++){
        if (j>0)
           comma=",";
        else
           comma="";
        s4+=comma+node.get(j).x;
        s5+=comma+node.get(j).y;
        s6+=comma+node.get(j).z;
    }
    System.out.println("X=["+s4+"];");
    System.out.println("Y=["+s5+"];");
    System.out.println("Z=["+s6+"];");
    System.out.println("plot3(X,Y,Z,'--rs','LineWidth',2,...");
    System.out.println("    'MarkerEdgeColor','k',...");
    System.out.println("    'MarkerFaceColor','g',...");
    System.out.println("    'MarkerSize',10)");
    System.out.println("axis square");
    System.out.println("grid on");
}
//////public static void main(String[] args){
//////generate12guide();
//////ArrayList<ArrayList> myArrays=new ArrayList<ArrayList>();
//////
//////myArrays.add(A);
//////myArrays.add(B);
//////
//////ArrayList<Node3d> A1=new ArrayList<Node3d>();
//////ArrayList<Node3d> B1=new ArrayList<Node3d>();
//////A1=myArrays.get(0);
//////B1=myArrays.get(1);
//////
//////printAarry(A);
//////printAarry(B);
//////
//////printAarry(A1);
//////printAarry(B1);
//////
////////printAarry(C);
////////printAarry(D);
////////printAarry(E);
////////printAarry(F);
////////printAarry(G);
////////printAarry(H);
////////printAarry(I);
////////printAarry(J);
////////printAarry(K);
////////printAarry(L);
////////////    generatePoints();
////////////    int counter=0;
////////////
////////////    for (int i=0;i<12;i++){
////////////    int[] X=new int[12];int[] Y=new int[12];int[] Z=new int[12];
////////////    String s4="",s5="",s6="";
////////////    String comma;
////////////    for(int j=0;j<12;j++){
////////////    X[j]=relativePoints.get(counter).x;
////////////    Y[j]=relativePoints.get(counter).y;
////////////    Z[j]=relativePoints.get(counter).z;
////////////        if (j>0)
////////////           comma=",";
////////////       else
////////////           comma="";
////////////        s4+=comma+X[j];
////////////        s5+=comma+Y[j];
////////////        s6+=comma+Z[j];
////////////    counter++;
////////////    }
////////////    System.out.println("X=["+s4+"];");
////////////    System.out.println("Y=["+s5+"];");
////////////    System.out.println("Z=["+s6+"];");
////////////    System.out.println("plot3(X,Y,Z,'--rs','LineWidth',2,...");
////////////    System.out.println("    'MarkerEdgeColor','k',...");
////////////    System.out.println("    'MarkerFaceColor','g',...");
////////////    System.out.println("    'MarkerSize',10)");
////////////    System.out.println("axis square");
////////////    System.out.println("grid on");
////////////    }
//////
//////}


}
