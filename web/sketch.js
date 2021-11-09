var api = new WSApi();

var x = 275;
var y = 275;
var t = 0;
var z = 1;
var r = 0;
var centerX;
var centerY;
var updateReady = true;
var speed = [0,0,0,0];

function setup() {
  createCanvas(windowWidth-50,windowHeight-50);
  centerX = (windowWidth-50)/2;
  centerY = (windowHeight-50)/2;
  x = centerX;
  y = centerY;
}

function drawBlade(x, y, direction, speed) {
  push();
  if (speed > 0) {
    colorMode(HSB);
    fill(color(210, 50*speed, 88));
  }
  translate(x,y);
  circle(0,0,20);
  rotate(PI / 3.0 + t*direction*speed);
  translate(-2.5, -10);
  rect(0,0, 5,20);
  pop();
}

function keyPressed() {
  console.log(key);
    api.sendCommand("keydown",{key: key, keyCode: key.charCodeAt(0)});
}

function keyReleased() {
    api.sendCommand("keyup",{key: key, keyCode: key.charCodeAt(0)});
}

function draw() {


	if (updateReady) {
		api.sendCommand("update", {delta: 0.0, simSpeed: 1.0}).then(function(updateData) {
			updateReady = true;
      x = centerX + updateData.pos[0]*100;
      y = centerY - updateData.pos[1]*100;
			z = updateData.pos[2] + 1;
      speed = updateData.propSpeed;
		});
		updateReady = false;
	}

  t = t + 0.1;

  background(220);
  
  translate(x, y);
  translate(16,25);
  scale(z);
  rotate(r);
  translate(-16,-25);
  fill(color(0, 0, 0));
  rect(11, 10, 10, 30);
  fill(color(255, 255, 255));

  
  drawBlade(0,0,-1,speed[0]);
  drawBlade(30,0,1,speed[1]);
  drawBlade(0,50,1,speed[2]);
  drawBlade(30,50,-1,speed[3]);
  drawVector();
}

function drawVector() {
  translate(x, y);

}
