package board;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.io.Serializable;
import java.util.ArrayList;

public class Stroke implements Serializable {
	Color color = Color.black;
	public ArrayList<Vector2D> points;

	Stroke() {
		points = new ArrayList<Vector2D>();
	}

	public void addPoint(Point p) {
		points.add(new Vector2D(p));
	}

	public void paint(Graphics g) {

		if (color == null)
			return;
		g.setColor(color);
		for (int i = 0; i < points.size() - 1; i++) {
			Vector2D p0 = (Vector2D) points.get(i);
			Vector2D p1 = (Vector2D) points.get(i + 1);
			if (p0 != null && p1 != null) {
				g.drawLine((int) p0.x, (int) p0.y, (int) p1.x, (int) p1.y);
			}
		}
		if (points.size() == 1) {
			Vector2D p0 = (Vector2D) points.get(0);
			g.drawLine((int) p0.x, (int) p0.y, (int) p0.x, (int) p0.y);
		}
	}

	public boolean intersects(Stroke stroke) {
		for (int i = 0; i < points.size() - 1; i++) {
			Vector2D p0 = (Vector2D) points.get(i);
			Vector2D p1 = (Vector2D) points.get(i + 1);
			if (stroke.intersects(p0, p1))
				return true;
		}
		return false;
	}

	public boolean intersects(Vector2D p0, Vector2D p1) {
		for (int i = 0; i < points.size() - 1; i++) {
			Vector2D q0 = (Vector2D) points.get(i);
			Vector2D q1 = (Vector2D) points.get(i + 1);
			if (intersects(p0.x, p0.y, p1.x, p1.y, q0.x, q0.y, q1.x, q1.y))
				return true;
		}
		return false;

	}

	// ü•ª‚Æ‚µ‚ÄŒð“_‚ ‚é‚©B
	public boolean intersects(
			double x1, double y1, double x2, double y2,
			double xx1, double yy1, double xx2, double yy2) {

		// a * x + b * y + c = 0
		double a0, b0, c0,
				a1, b1, c1;
		a0 = y1 - y2;
		b0 = x2 - x1;
		c0 = y2 * x1 - x2 * y1;
		a1 = yy1 - yy2;
		b1 = xx2 - xx1;
		c1 = yy2 * xx1 - xx2 * yy1;

		if (((a0 * xx1 + b0 * yy1 + c0) * (a0 * xx2 + b0 * yy2 + c0) <= 0) &&
				((a1 * x1 + b1 * y1 + c1) * (a1 * x2 + b1 * y2 + c1) <= 0))
			return true;
		else
			return false;
	}

}
