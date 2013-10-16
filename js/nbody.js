/* N-Body Simulation Usage
 * 
 * INITIALIZATION/RESET
 * var NumBodies = <insert number>;
 * var bodies = new Array(NumBodies);
 * for(var i = 0; i < NumBodies; i++)
 * 	bodies[i] = new GBody(new Vector3(<x0>,<y0>,<z0>),new Vector3(<vx0>,<vy0>,<vz0>),<mass>);
 * bodyField = new GBodyField(bodies);
 *
 * UPDATE CODE
 * bodyField.update(dt);
 *
 */

 var rk4 = function(a,r,v,dt){
        var a0    = a*r;
        var a1    = a*(r + 0.5*dt*v + 0.125*dt*dt*a0);
        var a2    = a*(r +     dt*v + 0.500*dt*dt*a1);
        var new_r =    r +     dt*v + ((a0+2*a1)*dt*dt)/6;
        return {r: new_r, v: v};
}

var euler = function(a,r,v,dt){
        var new_v = v+a*dt;
        var new_r = r+v*dt+0.5*a*dt*dt;
        return {r: new_r, v: new_v}
}

var heun = function(a,r,v,dt){
        //var new_v = v+2*a*dt;
        var new_r = r+2*v*dt;
        return {r: new_r, v: v};
}
 
 var numericalMethod = euler;
 
 var gravityDistanceThreshold = 0.1;
 var G = .00001;
 //var G = 0.0000000000667383
 
//Simple Vector Object
function Vector3(X,Y,Z){
	this.x = X;
	this.y = Y;
	this.z = Z;
}

//Body field controls all the gravitational bodies
function GBodyField(GBodies){
	this.bodies = GBodies;
	this.update = function(dt){ //Update given a dt (this is where the main computation for the simulation happens)
		var forceArray = new Array(this.bodies.length); //Array of computed forces
		for(var i = 0; i < this.bodies.length;i++){ //Compute each force
			forceArray[i] = new Array(this.bodies.length-i);
			var netForce = new Vector3(0,0,0); //Generate a net force for an object
			for(var j = 0; j < this.bodies.length;j++){ //For each object
				if(j<i){ //If we have computed this force before, look it up and rescale it to this current object mass
					netForce.x += (forceArray[j][i].x/(this.bodies[i].pos.x-this.bodies[j].pos.x))*(this.bodies[j].pos.x-this.bodies[i].pos.x);
					netForce.y += (forceArray[j][i].y/(this.bodies[i].pos.y-this.bodies[j].pos.y))*(this.bodies[j].pos.y-this.bodies[i].pos.y);
					netForce.z += (forceArray[j][i].z/(this.bodies[i].pos.z-this.bodies[j].pos.z))*(this.bodies[j].pos.z-this.bodies[i].pos.z);
				}
				else if (j>i){ //If we haven't computed this force before, add it to the force array
					forceArray[i][j] = this.bodies[i].calculateForce(this.bodies[j]);
					netForce.x += forceArray[i][j].x;
					netForce.y += forceArray[i][j].y;
					netForce.z += forceArray[i][j].z;
				}
			}
			//Runge-Kutta
			var new_x_vals = numericalMethod(netForce.x/this.bodies[i].m,this.bodies[i].pos.x,this.bodies[i].vel.x,dt);
			var new_y_vals = numericalMethod(netForce.y/this.bodies[i].m,this.bodies[i].pos.y,this.bodies[i].vel.y,dt);
			var new_z_vals = numericalMethod(netForce.z/this.bodies[i].m,this.bodies[i].pos.z,this.bodies[i].vel.z,dt);
			this.bodies[i].pos.x = new_x_vals.r;
			this.bodies[i].vel.x = new_x_vals.v;
			this.bodies[i].pos.y = new_y_vals.r;
			this.bodies[i].vel.y = new_y_vals.v;
			this.bodies[i].pos.z = new_z_vals.r;
			this.bodies[i].vel.z = new_z_vals.v;
		}
	}
}

//An object representing a single gravitational object
function GBody(Position, Velocity, Mass){
	this.pos = Position;
	this.vel = Velocity;
	this.m = Mass;
	this.calculateForce = function(body){ //Simple force calculation function
		//Calculate the changes in x and y direction and distance
		var dx = body.pos.x-this.pos.x;
		var dy = body.pos.y-this.pos.y;
		var dz = body.pos.z-this.pos.z;
		var d = Math.sqrt(dx*dx+dy*dy+dz*dz);
		//If the distances is too close, do not calculate force (prevents singularities)
		if(d<gravityDistanceThreshold) return new Vector3(0,0,0);
		//Create a force vector and return it
		var Fmodified = (this.m*body.m*G)/(d*d*d); //third 'd' is for unit vector
		dx*=Fmodified;
		dy*=Fmodified;
		dz*=Fmodified;
		return new Vector3(dx,dy,dz);
	}
}