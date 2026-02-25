import cairo

surface = cairo.SVGSurface("cairo_test2.svg", 500, 500)
ctx = cairo.Context(surface)
m = cairo.Matrix(yy=-1, y0=250, x0=250)
ctx.transform(m)
ctx.scale(30.0, 30.0)

ctx.save()
ctx.set_source_rgb(0.95, 0.95, 0.95)
ctx.paint()
ctx.restore()

ctx.new_path()
ctx.move_to(6.5, 0)
ctx.line_to(4.15, 0)
ctx.set_source_rgba(1.0, 0.0, 0.0, 1.0)
ctx.set_line_width(0.01) # test small line width
ctx.stroke()
surface.finish()
