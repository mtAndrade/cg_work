let ADD = 0.01;
let shape;
    
let createGeometry = function() {
    let geometry = new THREE.Geometry();

    geometry.vertices.push( new THREE.Vector3(3, 0, 0));
    geometry.vertices.push( new THREE.Vector3(0, 5, 0));
    geometry.vertices.push( new THREE.Vector3(0, 0, 2));
    geometry.vertices.push( new THREE.Vector3(1, 2, -2));
    
    geometry.faces.push(new THREE.Face3(0, 1, 2));
    geometry.faces.push(new THREE.Face3(1, 2, 3));
    geometry.faces.push(new THREE.Face3(0, 2, 3));
    
    
    let material = new THREE.MeshBasicMaterial(
                                    {color: 0xffffff,
                                    side:THREE.DoubleSide,     
                                    wireframe: true});
    shape = new THREE.Mesh( geometry, material );
    scene.add(shape);
};
