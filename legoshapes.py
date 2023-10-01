from svgelements import *;
import sys;
import bisect
import math

from types import MethodType

def split_cubic(curve, position):
    p00 = curve.start
    p01 = curve.control1
    p02 = curve.control2
    p03 = curve.end
    p10 = Point.towards(p00, p01, position)
    p11 = Point.towards(p01, p02, position)
    p12 = Point.towards(p02, p03, position)
    p20 = Point.towards(p10, p11, position)
    p21 = Point.towards(p12, p11, position)
    p30 = Point.towards(p20, p21, position)
    low = CubicBezier(p00, p10, p20, p30, smooth=False);
    high = CubicBezier(p30, p21, p12, p03, smooth=False);
    return low,high

def cubic_bezier_linear_error(curve):
    # Rotate and translate so start is at (0,0) and end is (x,0)
    end = curve.end - curve .start
    end_len = abs(end)
    s = -end.y / end_len
    c = end.x  / end_len
    m = Matrix(c,s,-s,c,
               -curve.start.x * c + curve.start.y * s,
               -curve.start.x * s - curve.start.y * c);
    flat = curve*m;
    xmin, ymin, xmax, ymax = flat.bbox()
    return max(-xmin, xmax - flat.end.x), max(ymax, -ymin)


def recursive_split(curve, start, end, max_error):
    error = cubic_bezier_linear_error(curve)
    if error[0] > 0 or error[1] > max_error:
        low, high = split_cubic(curve, 0.5)
        mid = (start + end) / 2;
        low_splits = recursive_split(low, start, mid, max_error)
        high_splits = recursive_split(high, mid, end, max_error)
        return low_splits + high_splits
    else:
        return [curve.start]

def cubic_polygon(self):
    polygon = recursive_split(self, 0,1, 1e-2)
    polygon.append(self.end)
    return polygon
    

CubicBezier.polygon = cubic_polygon;

def path_seg_polygon(self):
    if self.end == None or self.start == None or self.end == self.start:
        return [];
    else:
        return [self.start, self.end]
    
PathSegment.polygon = path_seg_polygon

def intersect_line_circle(line1, line2, center, radius):
    a0 = line1.x - center.x
    a1 = line2.x - line1.x
    b0 = line1.y - center.y
    b1 = line2.y - line1.y
    c = a1*a1 + b1*b1
    d = a0*b1 - a1*b0
    e = b0*b1 + a0*a1
    g = c*radius*radius - d*d
    if g <= 0: return None, None
    f = sqrt(g)
    s1 = (-e - f) / c
    s2 = (-e + f) / c
    p3_1 = Point.towards(line1, line2, s1) if s1 >= 0 and s1 < 1 else None
    p3_2 = Point.towards(line1, line2, s2) if s2 >= 0 and s2 < 1 else None
    return p3_1, p3_2
    
"""p rotated 90 degrees"""
def rot90(p):
    return Point(-p.y, p.x)

"""Change length of p to 1"""
def unit(p):
    return p * (1 / abs(p))

def poly_offset(poly, offset):
    p_n1 = poly[0]
    p_n2 = None
    offseted = []
    for point in poly[1:]:
        dir = unit(rot90(point - p_n1))
        if p_n2 != None:
            dir = unit(unit(rot90(p_n1 - p_n2)) + dir)
        offseted.append(dir * offset + p_n1)
        p_n2 = p_n1;
        p_n1 = point
    dir = unit(rot90(p_n1 - p_n2))
    offseted.append(dir * offset + p_n1)
    return offseted

def poly_bbox(poly):
    xmin=poly[0].x
    xmax=poly[0].x
    ymin=poly[0].y
    ymax=poly[0].y
    
    for p in poly:
        xmax = max(xmax, p.x)
        xmin = min(xmin, p.x)
        ymax = max(ymax, p.y)
        ymin = min(ymin, p.y)
    return xmin, ymin, xmax, ymax

def append_seg(segs,seg, seg_inc):
    if len(seg) >= 2:
        if not seg_inc:
            seg.reverse()
        segs.append(seg)
            
"""Split the polygon into segments with increasing y coordinates. The
polygon is assumed to be closed, i.e. the first coordinate must be equal to the last"""
def split_on_direction(poly):
    seg_inc = poly[0].y < poly[1].y # The current segment has increasing y
    prev_point = poly[0]
    seg = [poly[0]]
    segs = []
    for point in poly[1:]:
        if prev_point.y != point.y:
            inc = prev_point.y < point.y
            if inc != seg_inc:
                append_seg(segs, seg, seg_inc)
                seg=[prev_point]
                seg_inc = inc
        else:
            append_seg(segs, seg, seg_inc)
            seg = []
        seg.append(point)
        prev_point = point

    append_seg(segs, seg, seg_inc)

    return segs

def point_in_poly(segs, point):
    count = 0
    for seg in segs:
        i = bisect.bisect_left(seg, point.y, key = lambda a: a.y)
        if i > 0 and i < len(seg):
            x2 = seg[i-1].x
            y2 = seg[i-1].y
            x1 = seg[i].x
            y1 = seg[i].y

            if (x1 -point.x) * (y1 - y2) + (x2-x1)*(y1-point.y) > 0:
                count += 1
    return (count & 1) == 1

def remove_short(poly, min):
    prev = None
    for p in poly:
        if prev == None or abs(prev - p) >= min:
            prev = p
            yield p


if len(sys.argv) <= 2:
    out = sys.stdout
else:
    out = open(sys.argv[2], "w")

WALL_THICKNESS=1.2
TOP_THICKNESS=1
STUD_DIST=8
ANTISTUD_CLEARANCE=3.5+WALL_THICKNESS
BRICK_HEIGHT=3.2
elems = SVG.parse(sys.argv[1]).elements();
for elem in elems:
    path = None
    if isinstance(elem, Path):
        if len(elem) != 0:
            path = Path(elem)
    elif isinstance(elem, Shape):
        path = Path(elem)
        path.reify()  # In some cases the shape could not have reified, the path must.

    if path != None:
        print(path)
        path_poly = []
        for seg in path.segments():
            poly = seg.polygon()
            path_poly.extend(poly);

        path_poly = list(remove_short(path_poly, 1e-6))
        if len(path_poly) >= 2:
            offseted = poly_offset(path_poly,5)                
        segs = split_on_direction(path_poly)
        out.write("include <shape_parts.scad>\n"  )
        out.write("module shape() {\n");
        out.write("polygon(points=["+",".join("[%f, %f]" % (p.x,p.y) for p in path_poly)+ "]);\n"  )
        out.write("}\n"  )
        
        out.write("translate([0,0,%f]) {\n"%(BRICK_HEIGHT-TOP_THICKNESS))
        out.write("linear_extrude(height=%f) {\n"%TOP_THICKNESS  )
        out.write("shape();\n"  )
        out.write("}\n"  )
        out.write("}\n"  )
        
        out.write("difference() {\n"  )
        out.write("linear_extrude(height=%f) {\n"%BRICK_HEIGHT  )
        out.write("difference() {\n"  )
        out.write("shape();\n"  )
        out.write("offset(delta="+str(-WALL_THICKNESS)+") {shape();}\n"  )
        out.write("}\n"  )
        out.write("}\n"  )
        
        xmin, ymin, xmax, ymax = poly_bbox(path_poly)
        xmin = math.floor((xmin-STUD_DIST/2) / STUD_DIST) * STUD_DIST
        ymin = math.floor((ymin-STUD_DIST/2) / STUD_DIST) * STUD_DIST
        y = ymin
        while y < ymax:
            x = xmin
            while x < xmax:
                 out.write("translate([%f,%f,0]) {negative_stud();}\n"
                           %(x+STUD_DIST/2,y+STUD_DIST/2))
                 x += STUD_DIST
            y += STUD_DIST
        out.write("}\n"  )

        y = ymin
        while y < ymax:
            x = xmin
            while x < xmax:
                if point_in_poly(segs, Point(x,y)):
                    prev_p = path_poly[0]
                    intersected = False
                    for p in path_poly[1:]:
                        if intersect_line_circle(prev_p, p, Point(x,y), ANTISTUD_CLEARANCE) != (None,None):
                            intersected = True
                            break
                        prev_p = p
                    if not intersected:
                         out.write("translate([%f,%f]) {antistud(%f);}\n"%(x,y, BRICK_HEIGHT))
                x += STUD_DIST
            y += STUD_DIST
        
