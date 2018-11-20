import svgwrite

start    = (100,100)
end      = (200,100)
radius   = 50

dwg      = svgwrite.Drawing('arc_test.svg')
arc_path = dwg.path(d='M %f , %f' % start, fill='none', stroke='black')
arc_path.push_arc(end, 0, (50,1), False, '+', True)

dwg.add(arc_path)

dwg.save()




