package board;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Panel;
import java.awt.Point;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;

class MainPanel extends Panel implements ButtonsInterface {

	MainPanel() {
		MouseDispatcher mouseDispatcher = new MouseDispatcher();
		addMouseListener(mouseDispatcher);
		addMouseMotionListener(mouseDispatcher);

	}

	public void initializeButtons(Buttons buttons) {
		buttons.addButton("quit");
		buttons.addButton("clear");
	}

	public void buttonPressed(String label) {
		if (label.equals("quit"))
			System.exit(0);
		if (label.equals("clear"))
			clear();
	}

	Image buffered_image = null;

	public void update(Graphics g) {
		paint(g);
	}

	public boolean repaint_all = true;

	public void repaint_all() {
		repaint_all = true;
		repaint();
	}

	public void paint(Graphics g) {

		// double buffering
		if (buffered_image == null)
			buffered_image = createImage(getSize().width, getSize().height);
		Graphics bg = buffered_image.getGraphics();

		// clear
		bg.setColor(Color.white);
		bg.fillRect(0, 0, getSize().width, getSize().height);

		// paint code!
		for (Oval oval : ovals)
			oval.paint(bg);
		// double buffering
		bg.dispose();

		g.drawImage(buffered_image, 0, 0, this);

		repaint_all = false;

	}

	ArrayList<Oval> ovals = new ArrayList<Oval>();

	//
	// event handling
	//
	Oval oval;

	public void start_stroke(Point p) {
		oval = new Oval(p);
		ovals.add(oval);
	}

	public void extend_stroke(Point p) {
		oval.extend(p);
		repaint();
	}

	public void finish_stroke(Point p) {
	}

	public void clear() {
		ovals.clear();
		repaint_all();
	}

	public int operation_status;
	public static final int NONE = 0;
	public static final int DRAWING = 1;
	public static final int RESIZING = 2;

	//
	// central event dispatcher
	//
	public class MouseDispatcher
			extends MouseAdapter implements MouseMotionListener {

		public void mouseMoved(MouseEvent e) {
		}

		public void mousePressed(MouseEvent e) {
			Point p = e.getPoint();
			boolean right_button = (e.getModifiers() & e.BUTTON3_MASK) != 0;

			operation_status = DRAWING;
			start_stroke(p);
		}

		public void mouseDragged(MouseEvent e) {
			Point p = e.getPoint();

			switch (operation_status) {
				case DRAWING:
					extend_stroke(p);
					break;
			}

		}

		public void mouseReleased(MouseEvent e) {
			Point p = e.getPoint();

			switch (operation_status) {
				case DRAWING:
					finish_stroke(p);
					break;
			}
		}

		public void mouseClicked(MouseEvent e) {

		}
	}

	public static Frame frame;

	public static void main(String args[]) {
		frame = new Frame("Main");
		frame.addWindowListener(
				new WindowAdapter() {
					public void windowClosing(WindowEvent e) {
						System.exit(0);
					}
				});
		frame.addKeyListener(
				new KeyAdapter() {
					public void keyPress(KeyEvent e) {
						if (e.getKeyCode() == e.VK_C)
							System.exit(0);
					}
				});

		Panel panel = new Panel();
		panel.setLayout(new BorderLayout());
		MainPanel mainPanel = new MainPanel();
		panel.add("Center", mainPanel);
		panel.add("South", new Buttons(mainPanel));

		frame.add("Center", panel);
		frame.setSize(600, 600);
		frame.setVisible(true);

	}
}
