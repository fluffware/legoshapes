$fn=100;
STUD_DIST=8;
STUD_HEIGHT=1.8;
NEG_STUD_DIAM=4.8;
ANTISTUD_OUTER_DIAM = 6.51;
ANTISTUD_INNER_DIAM = 4.8;
module antistud(height) {
  difference() {
    cylinder(d=ANTISTUD_OUTER_DIAM, h=height);
    translate([0,0,-1]) {
      cylinder(d=ANTISTUD_INNER_DIAM, h=height+2);
    }
  }
}

module negative_stud() {
  cylinder(d=NEG_STUD_DIAM, h=STUD_HEIGHT);
}
