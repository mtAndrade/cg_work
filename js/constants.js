
const M_PI =  3.14159265358979323846;
const TWOPI =  (2*M_PI);
const LEN = (x,y,z) =>{  return (Math.sqrt( (x)*(x) + (y)*(y) +(z)*(z) ) ) }
const clamp = (x) => {return  (Math.min( Math.max((x),0),1))}

const makeQuad = function(geometry, position, addFace, verts) {  
  geometry.vertices.push(position);
  if (addFace) {
    var index1 = geometry.vertices.length - 1;
    var index2 = index1 - 1;
    var index3 = index1 - verts;
    var index4 = index1 - verts - 1;
    
    geometry.faces.push(new THREE.Face3(index2, index3, index1));
    geometry.faces.push(new THREE.Face3(index2, index4, index3));
  }
};

let setBoxAndNormal = function(geometry, mesh, invert_orientation=false){
  if(showNormals){
    if(invert_orientation){
      for ( var i = 0; i < geometry.faces.length; i ++ ) {
          var face = geometry.faces[ i ];
          var temp = face.a;
          face.a = face.c;
          face.c = temp;
      }
      geometry.computeFaceNormals();
      geometry.computeVertexNormals();
      var faceVertexUvs = geometry.faceVertexUvs[ 0 ];
      for ( var i = 0; i < faceVertexUvs.length; i ++ ) {
          var temp = faceVertexUvs[ i ][ 0 ];
          faceVertexUvs[ i ][ 0 ] = faceVertexUvs[ i ][ 2 ];
          faceVertexUvs[ i ][ 2 ] = temp;
      }
    }
    
    geometry.computeFaceNormals()
    normals = new THREE.FaceNormalsHelper( mesh, 0.2, 0x00ff00, 1 );
    scene.add(normals)
  }
  if(showBox){
      geometry.computeBoundingBox()
      box = new THREE.BoxHelper(mesh, 0x000000 );
      scene.add( box );
  }
  if(axis){
    var axesHelper = new THREE.AxesHelper( 3 );
    scene.add( axesHelper );
  }
}