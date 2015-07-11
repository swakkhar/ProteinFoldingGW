/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

/**
 *
 * @author mar
 */
public class Node3d {

    public int x;
    public int y;
    public int z;
    Node3d()
    {
        this.x=0;
        this.y=0;
        this.z=0;
    }
    Node3d(int XCoordinate,int YCoordinate,int ZCoordinate)
    {
        this.x=XCoordinate;
        this.y=YCoordinate;
        this.z=ZCoordinate;
    }
    Node3d(int[] point)
    {
        this.x=point[0];
        this.y=point[1];
        this.z=point[2];
    }
    Node3d(Node3d node)
    {
        this.x=node.x;
        this.y=node.y;
        this.z=node.z;
    }
    public boolean isEqual(Node3d node){
    if (node.x==this.x && node.y==this.y && node.z==this.z)
        return true;
    else
        return false;
    }
}
