import svgwrite

dwg = svgwrite.Drawing('test.svg')
#dwg.add(dwg.line((0, 0), (0, 0), stroke=svgwrite.rgb(255, 0, 0, '%')))
dwg.add(dwg.circle(center=(1, 1), r=2))
#dwg.add(dwg.text('Test', insert=(5, 12), fill='red'))
dwg.save()

