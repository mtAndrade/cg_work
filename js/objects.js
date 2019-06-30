
let computeUVs = function(geometry) {
    if (!geometry.boundingBox) geometry.computeBoundingBox();

    let sz      = geometry.boundingBox.getSize(new THREE.Vector3());
    let center  = geometry.boundingBox.getCenter(new THREE.Vector3())
    let min     = geometry.boundingBox.min;

    if (geometry.faceVertexUvs[0].length == 0) {
        for (let i = 0; i < geometry.faces.length; i++) {
            geometry.faceVertexUvs[0].push([new THREE.Vector2(), new THREE.Vector2(), new THREE.Vector2()]);
        }
    }

    for (let i = 0; i < geometry.faces.length; i++) {
        let face    = geometry.faces[i];
        let faceUVs = geometry.faceVertexUvs[0][i]
        
        let va  = geometry.vertices[geometry.faces[i].a]
        let vb  = geometry.vertices[geometry.faces[i].b]
        let vc  = geometry.vertices[geometry.faces[i].c]
        let vab = new THREE.Vector3().copy(vb).sub(va)
        let vac = new THREE.Vector3().copy(vc).sub(va)
        
        let vcross = new THREE.Vector3().copy(vab).cross(vac);
        
        vcross.set(Math.abs(vcross.x), Math.abs(vcross.y), Math.abs(vcross.z))
        
        let majorAxis = vcross.x > vcross.y ? (vcross.x > vcross.z ? 'x' : vcross.y > vcross.z ? 'y' : vcross.y > vcross.z) : vcross.y > vcross.z ? 'y' : 'z'
        
        let uAxis = majorAxis == 'x' ? 'y' : majorAxis == 'y' ? 'x' : 'x';
        let vAxis = majorAxis == 'x' ? 'z' : majorAxis == 'y' ? 'z' : 'y';
        faceUVs[0].set((va[uAxis] - min[uAxis]) / sz[uAxis], (va[vAxis] - min[vAxis]) / sz[vAxis])
        faceUVs[1].set((vb[uAxis] - min[uAxis]) / sz[uAxis], (vb[vAxis] - min[vAxis]) / sz[vAxis])
        faceUVs[2].set((vc[uAxis] - min[uAxis]) / sz[uAxis], (vc[vAxis] - min[vAxis]) / sz[vAxis])
    }
    
    geometry.elementsNeedUpdate = geometry.verticesNeedUpdate = true;
}

let createCone = function(h, r, vs, rs){
    const hr = h / r;
    
    const geometry = new THREE.Geometry()
    let points = []
    let norms = []

    for (i = 0; i < vs; i++) {
       for (j = 0; j <= rs; j++) {
          for (k = 1; k >= 0; k--) {
             z = h * (i+k) / vs;              
             u = (h-z) / hr; // cone is upside down
             
             th = j * TWOPI / rs;
             x = u * Math.cos(th);
             y = u * Math.sin(th);
 
             nx = x * hr; ny = y * hr; nz = u > 0 ? u : 1; // cone is upside down
             ln = LEN(nx,ny,nz);
             nx /= ln; ny /= ln; nz /= ln;
             th = clamp(th/TWOPI); 

             const position = new THREE.Vector3(x, y, z);
             const addFace = (i >= 0 ) && (j > 0) && (k==0);
             makeQuad(geometry, position, addFace,   2);             
          }
       }
    }
    computeUVs(geometry);
    
    cone = new THREE.Mesh( geometry, defaultMaterial );
    
    setBoxAndNormal(geometry, cone, invert_faces)

    scene.add(cone);

    return cone
}

let createCylinder = function(h, r, vs, rs) { 
    let geometry = new THREE.Geometry();    
    let points = []

    for (i = 0; i < vs; i++) {       // generate stacks
       z0 = h * i / vs;              // stack bottom z
       z1 = h * (i+1) / vs;          // stack top z
       
       for (j = 0; j <= rs; j++) {   // generate slices
            th = j * TWOPI / rs;     // slice left theta
            x = r * Math.cos(th);
            y = r * Math.sin(th);
 
            ln = LEN(x,y,0);
            nx = x/ln; ny = y/ln;
            th = clamp(th/TWOPI);
 
            points.push(x,y,z1)
            points.push(x,y,z0)
            const position = new THREE.Vector3(x, y, z1);
            const position2 = new THREE.Vector3(x, y, z0);
            const addFace = (i > 0 ) && (j > 0) ;
            makeQuad(geometry, position, addFace,   2);
            makeQuad(geometry, position2, addFace,   2);
       }
    }
    computeUVs(geometry);
    
    cylinder = new THREE.Mesh( geometry, defaultMaterial );
    
    setBoxAndNormal(geometry, cylinder, invert_faces)

    scene.add(cylinder);
    return cylinder
}

let createParaboloid = function(h, r, vs, rs) {
    let geometry = new THREE.Geometry();    
    let points = []

    const hr2 = h/(r*r);

    for (i = 0; i < vs; i++) {
       for (j = 0; j <= rs; j++) {
          for (k = 1; k >= 0; k--) {
             z = h * (i+k) / vs;              
             u = Math.sqrt((h-z)/hr2); // paraboloid is upside down
 
             th = j * TWOPI / rs;
             x = u * Math.cos(th);
             y = u * Math.sin(th);
 
             nx = 2 * x * hr2; ny = 2 * y * hr2; nz = 1; // paraboloid is upside down
             ln = LEN(nx,ny,nz);
             nx /= ln; ny /= ln; nz /= ln;
             th = clamp(th/TWOPI);
 
             points.push(x, y, z)
             const position = new THREE.Vector3(x, y, z);
             const addFace = (i >= 0 ) && (j > 0) && (k==0);
             makeQuad(geometry, position, addFace,   2);
          }
       }
    }
    computeUVs(geometry);
    
    paraboloid = new THREE.Mesh( geometry, defaultMaterial );
    
    setBoxAndNormal(geometry, paraboloid, invert_faces)
    scene.add(paraboloid);
    return paraboloid    
 }

let createHyperboloid = function(h, a, r, vs, rs) {
    let geometry = new THREE.Geometry();
    let points = []
    
    const c = Math.sqrt( h*h / (4 * (r*r/(a*a) - 1)) );
    if ( r < a ) {
       console.log( "Hyperboloid: r < a -> invalid!\n" );
       return;
    }
    
    for (i = -vs; i < vs; i++) {
        for (j = 0; j <= rs; j++) {
            for (k = 1; k >= 0; k--) {
                z = h * (i+k) / vs;              
                u = z/c;
                u = a * Math.sqrt(1+u*u);

                th = j * TWOPI / rs;
                x = u * Math.cos(th);
                y = u * Math.sin(th);

                nx = x * c; ny = y * c; nz = -z * a * a / c;
                ln = LEN(nx,ny,nz);
                nx /= ln; ny /= ln; nz /= ln;
                th = clamp(th/TWOPI);

                points.push(x, y, z);
                const position = new THREE.Vector3(x, y, z);
                const addFace = (i >= -vs ) && (j > 0) && (k==0);
                makeQuad(geometry, position, addFace,   2);
            }
        }
    }
    computeUVs(geometry);

    hyperboloid = new THREE.Mesh( geometry, defaultMaterial );
    
    setBoxAndNormal(geometry, hyperboloid, invert_faces)
    scene.add(hyperboloid);
    return hyperboloid
}
 
let createBufferHyperboloid = function(h, a, r, vs, rs) {
    let buffer_geometry = new THREE.BufferGeometry();
    let points = []
    
    const c = Math.sqrt( h*h / (4 * (r*r/(a*a) - 1)) );
    if ( r < a ) {
       console.log( "Hyperboloid: r < a -> invalid!\n" );
       return;
    }
    
    for (i = -vs; i < vs; i++) {
        for (j = 0; j <= rs; j++) {
            for (k = 1; k >= 0; k--) {
                z = h * (i+k) / vs;              
                u = z/c;
                u = a * Math.sqrt(1+u*u);

                th = j * TWOPI / rs;
                x = u * Math.cos(th);
                y = u * Math.sin(th);

                nx = x * c; ny = y * c; nz = -z * a * a / c;
                ln = LEN(nx,ny,nz);
                nx /= ln; ny /= ln; nz /= ln;
                th = clamp(th/TWOPI);

                points.push(x, y, z);
            }
        }
    }
    const vertices = new Float32Array( points )
    buffer_geometry.addAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );

    hyperboloid = new THREE.Mesh( buffer_geometry, defaultMaterial );
    hyperboloid.drawMode = THREE.TriangleStripDrawMode;

    if(normals){
        buffer_geometry.computeFaceNormals()
        var helper = new THREE.FaceNormalsHelper( hyperboloid, 0.2, 0x00ff00, 1 );
        scene.add(helper)
    }
    
    scene.add(hyperboloid);
}

let createHyperbolicParaboloid = function(w, h, ws, hs) {
    let geometry = new THREE.Geometry();    
    let points = []
    
    for (i = -hs; i < hs; i++) {
       for (j = -ws; j <= ws; j++) {
          for (k = 1; k >= 0; k--) {
             u = j*w/ws;
             v = (i+k)*h/hs;	 
 
             x = u;
             y = v;
             z = u*v;
 
             nx = y; ny = x; nz = -1; 
             ln = LEN(nx,ny,nz);
             nx /= ln; ny /= ln; nz /= ln;

             points.push(x, y, z)
             const position = new THREE.Vector3(x, y, z);
             const addFace = (i >=-hs ) && (j > -ws) && (k==0);
             makeQuad(geometry, position, addFace,   2);
          }
       }
    }
    computeUVs(geometry);
    
    hyperbolic = new THREE.Mesh( geometry, defaultMaterial );
    
    setBoxAndNormal(geometry, hyperbolic, invert_faces)
    scene.add(hyperbolic);
    return hyperbolic
 }

let createSphere = function( r1, r2, r3, numc, numt) {
    let geometry = new THREE.Geometry();    
    for (i = 0; i < numc; i++) {
       for (j = 0; j <= numt; j++) {
          for (k = 1; k >= 0; k--) {
             s = (i + k); 
 
             u = j*TWOPI/numt;
             v = s*M_PI/numc;
 
             cu = Math.cos(u);
             cv = Math.cos(v);
             su = Math.sin(u);
             sv = Math.sin(v);
 
             // (x,y,z) = f(u,v)
             x = r1 * cu * sv;
             y = r2 * su * sv;
             z = r3 * cv;
 
             // df/du x df/dv
             nx = r2*r3*r2*r3 * x;
             ny = r1*r3*r1*r3 * y;
             nz = r1*r2*r1*r2 * z;
 
             ln = LEN(nx,ny,nz);
             nx /= ln; ny /= ln; nz /= ln;
 
             u = clamp(u/TWOPI);
             v = clamp(v/M_PI); 

             const position = new THREE.Vector3(x, y, z);
             const addFace = (i >= 0 ) && (j > 0) && (k==0);
             makeQuad(geometry, position, addFace,   2);
          }
       }
    }
    computeUVs(geometry);
    
    sphere = new THREE.Mesh( geometry, defaultMaterial );

    setBoxAndNormal(geometry, sphere, invert_faces)
    scene.add(sphere)
    return sphere
 }

let createTorus = function(r1, r2, numc, numt){
    let points = []
    let geometry = new THREE.Geometry()

    for (i = 0; i < numc; i++) {
        for (j = 0; j <= numt; j++) {
            for (k = 1; k >= 0; k--) {
                s = (i + k); 
                t = j;
                u = t*TWOPI/numt;
                v = s*TWOPI/numc;

                cu = Math.cos(u);
                cv = Math.cos(v);
                su = Math.sin(u);
                sv = Math.sin(v);

                // (x,y,z) = f(u,v)
                x = (r1 + r2 * cv) * cu;
                y = (r1 + r2 * cv) * su;
                z = r2 * sv;

                // df/du x df/dv
                nx = cv * cu;
                ny = cv * su;
                nz = sv;

                ln = LEN(nx,ny,nz);
                nx /= ln; ny /= ln; nz /= ln;

                u = clamp(u/TWOPI);
                v = clamp(v/TWOPI);

                points.push(x, y, z);

                const position = new THREE.Vector3(x, y, z);
                const addFace = (i >= 0 ) && (j > 0) && (k==0);
                makeQuad(geometry, position, addFace,   2);
            }
        }
    }
    computeUVs(geometry);
    
    torus = new THREE.Mesh( geometry, defaultMaterial );

    setBoxAndNormal(geometry, torus, invert_faces)
    scene.add(torus)
    return torus
}

 
