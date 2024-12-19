package board;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;

public class Oval {
	Point center;
	int w;
	int h;
	Point start_center;
	int start_w;
	int start_h;
	boolean selected = false;
	double angle = 0; // 回転角度
	double start_angle = 0;

	public Oval(Point p) {
		center = p;
	}

	public Oval(Oval oval) {
		center = new Point(oval.center);
		w = oval.w;
		h = oval.h;
		angle = oval.angle;
	}

	public void extend(Point p) {
		w = Math.abs(p.x - center.x) * 2;
		h = Math.abs(p.y - center.y) * 2;
	}

	public void start(Point p) {
		start_center = new Point(center);
		start_w = w;
		start_h = h;
		start_angle = angle;
	}

	public void resize(Point p, Point start_point, int selected_handle) {
		Point[] handles = getHandles();
		double cos = Math.cos(angle);
		double sin = Math.sin(angle);
		int dx, dy, dw, dh;
		dx = p.x - start_point.x;
		dy = p.y - start_point.y;
		dw = (int) (dx * cos + dy * sin);
		dh = (int) (-dx * sin + dy * cos);
		switch (selected_handle) {
			case 0: // w方向のみ変更
				center.x = (int) (start_center.x + dw / 2 * cos);
				center.y = (int) (start_center.y + dw / 2 * sin);
				w = start_w - dw;
				break;
			case 1: // w方向のみ変更
				center.x = (int) (start_center.x + dw / 2 * cos);
				center.y = (int) (start_center.y + dw / 2 * sin);
				w = start_w + dw;
				break;
			case 2: // h方向のみ変更
				center.x = (int) (start_center.x + dh / 2 * -sin);
				center.y = (int) (start_center.y + dh / 2 * cos);
				h = start_h - dh;
				break;
			case 3: // h方向のみ変更
				center.x = (int) (start_center.x + dh / 2 * -sin);
				center.y = (int) (start_center.y + dh / 2 * cos);
				h = start_h + dh;
				break;
			case 4: // w方向、h方向とも変更
				center.x = (int) (start_center.x + (dw / 2) * cos + (dh / 2) * -sin);
				center.y = (int) (start_center.y + (dw / 2) * sin + (dh / 2) * cos);
				w = start_w - dw;
				h = start_h - dh;
				break;
			case 5: // w方向、h方向とも変更
				center.x = (int) (start_center.x + (dw / 2) * cos + (dh / 2) * -sin);
				center.y = (int) (start_center.y + (dw / 2) * sin + (dh / 2) * cos);
				w = start_w + dw;
				h = start_h - dh;
				break;
			case 6: // w方向、h方向とも変更
				center.x = (int) (start_center.x + (dw / 2) * cos + (dh / 2) * -sin);
				center.y = (int) (start_center.y + (dw / 2) * sin + (dh / 2) * cos);
				w = start_w + dw;
				h = start_h + dh;
				break;
			case 7: // w方向、h方向とも変更
				center.x = (int) (start_center.x + (dw / 2) * cos + (dh / 2) * -sin);
				center.y = (int) (start_center.y + (dw / 2) * sin + (dh / 2) * cos);
				w = start_w - dw;
				h = start_h + dh;
				break;
		}
	}

	public void move(Point p, Point start_point) {
		center = new Point(start_center.x + p.x - start_point.x, start_center.y + p.y - start_point.y);
		// System.out.println("move");
	}

	public void rotate(Point p, Point start_point) {
		double dx = p.x - center.x;
		double dy = p.y - center.y;
		double angle1 = Math.atan2(dy, dx);
		dx = start_point.x - center.x;
		dy = start_point.y - center.y;
		double angle2 = Math.atan2(dy, dx);
		angle = start_angle + angle1 - angle2;
		// System.out.println("rotate");
	}

	public void paint(Graphics g) {
		Graphics2D g2d = (Graphics2D) g;
		AffineTransform origin = g2d.getTransform();
		if (selected) {
			g2d.setColor(Color.red);
		} else {
			g2d.setColor(Color.black);
		}

		// 楕円を定義
		Ellipse2D ellipse = new Ellipse2D.Double(center.x - w / 2, center.y - h / 2, w, h);

		// 回転用のAffineTransformを作成
		// originをコピーして、回転を加える
		AffineTransform at = new AffineTransform(origin);
		at.rotate(angle, center.x, center.y);

		// Graphics2Dに変換を設定
		g2d.setTransform(at);

		// 楕円を描画
		g2d.fill(ellipse);

		// 変換を元に戻す
		g2d.setTransform(origin);
		if (selected) {
			g2d.setColor(Color.gray);
			Point[] handles = getHandles();
			// 長方形で囲む
			g2d.drawLine(handles[4].x, handles[4].y, handles[5].x, handles[5].y);
			g2d.drawLine(handles[5].x, handles[5].y, handles[6].x, handles[6].y);
			g2d.drawLine(handles[6].x, handles[6].y, handles[7].x, handles[7].y);
			g2d.drawLine(handles[7].x, handles[7].y, handles[4].x, handles[4].y);
			for (int i = 0; i < 8; i++) {
				g2d.setColor(Color.white);
				g2d.fillOval(handles[i].x - 3, handles[i].y - 3, 6, 6);
				g2d.setColor(Color.black);
				g2d.drawOval(handles[i].x - 3, handles[i].y - 3, 6, 6);
			}
		}
	}

	public boolean contains(Point p) {
		double dx = p.x - center.x;
		double dy = p.y - center.y;
		double a = Math.cos(angle);
		double b = Math.sin(angle);
		double x = a * dx + b * dy;
		double y = -b * dx + a * dy;
		return x * x / (w * w / 4) + y * y / (h * h / 4) <= 1;
	}

	public Point[] getHandles() {
		Point[] handles = new Point[8];
		handles[0] = new Point(-w / 2, 0);
		handles[1] = new Point(w / 2, 0);
		handles[2] = new Point(0, -h / 2);
		handles[3] = new Point(0, h / 2);
		handles[4] = new Point(-w / 2, -h / 2);
		handles[5] = new Point(w / 2, -h / 2);
		handles[6] = new Point(w / 2, h / 2);
		handles[7] = new Point(-w / 2, h / 2);
		// rotate
		for (int i = 0; i < 8; i++) {
			double dx = handles[i].x;
			double dy = handles[i].y;
			double a = Math.cos(angle);
			double b = Math.sin(angle);
			handles[i].x = (int) (a * dx - b * dy) + center.x;
			handles[i].y = (int) (b * dx + a * dy) + center.y;
		}
		return handles;
	}

	public int selected_handle(Point p) {
		Point[] handles = getHandles();
		int area = 10;
		for (int i = 0; i < 8; i++) {
			if (Math.pow(handles[i].x - p.x, 2) + Math.pow(handles[i].y - p.y, 2) < Math.pow(area, 2)) {
				return i;
			}
		}
		return -1;
	}

	public void start_selected() {
		selected = true;
	}

	public void end_selected() {
		selected = false;
	}
}
