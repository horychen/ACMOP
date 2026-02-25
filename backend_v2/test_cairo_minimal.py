import cairo

# Create SVG
surface = cairo.SVGSurface("cairo_test.svg", 500, 500)
ctx = cairo.Context(surface)
m = cairo.Matrix(yy=-1, y0=250, x0=250)
ctx.transform(m)
ctx.scale(30.0, 30.0)

# Background
ctx.save()
ctx.set_source_rgb(0.95, 0.95, 0.95)
ctx.paint()
ctx.restore()

# Draw line
ctx.new_path()
ctx.set_source_rgba(0,0,0,1.0)
ctx.set_line_width(1.0/30.0)
ctx.move_to(0,0)
ctx.line_to(10,10)
ctx.stroke()

surface.finish()
print("Done")
