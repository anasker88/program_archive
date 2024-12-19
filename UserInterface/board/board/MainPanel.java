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

import javax.swing.JFrame;
import javax.swing.JOptionPane;

class MainPanel extends Panel implements ButtonsInterface {

	MainPanel() {
		MouseDispatcher mouseDispatcher = new MouseDispatcher();
		addMouseListener(mouseDispatcher);
		addMouseMotionListener(mouseDispatcher);

	}

	public void initializeButtons(Buttons buttons) {
		buttons.addButton("quit");
		buttons.addButton("clear");
		buttons.addButton("delete");
		buttons.addButton("copy");
	}

	public void buttonPressed(String label) {
		if (label.equals("quit"))
			System.exit(0);
		if (label.equals("clear"))
			clear();
		if (label.equals("delete"))
			delete();
		if (label.equals("copy"))
			copy();
	}

	Image buffered_image = null;

	public void update(Graphics g) {
		paint(g);
	}

	public boolean repaint_all = true;

	public void repaint_all() {
		repaint_all = true;
		repaint();
		// System.out.println(operation_status);
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

	public void delete() {
		if (operation_status == SELECTING) {
			if (selected_oval != null) {
				ovals.remove(selected_oval);
				selected_oval = null;
				repaint_all();
			}
		} else {
			JFrame frame = new JFrame();
			JOptionPane.showMessageDialog(frame, "Select an oval to delete.");
		}
		operation_status = NONE;
	}

	public void copy() {
		if (operation_status == SELECTING) {
			if (selected_oval != null) {
				Oval oval = new Oval(selected_oval);
				ovals.add(oval);
				oval.center.x += (int) (oval.w / 2);
				selected_oval.end_selected();
				selected_oval = oval;
				selected_oval.start_selected();
				repaint_all();
			}
		} else {
			JFrame frame = new JFrame();
			JOptionPane.showMessageDialog(frame, "Select an oval to copy.");
		}
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

	public Point start_point = null;

	public void start_resize(Point p) {
		start_point = p;
		selected_oval.start(p);
	}

	public void extend_resize(Point p) {
		selected_oval.resize(p, start_point, selected_handle);
		repaint_all();
	}

	public void finish_resize(Point p) {
		start_point = null;
	}

	public void start_move(Point p) {
		start_point = p;
		selected_oval.start(p);
	}

	public void extend_move(Point p) {
		selected_oval.move(p, start_point);
		repaint_all();
	}

	public void finish_move(Point p) {
		start_point = null;
	}

	public void start_rotate(Point p) {
		start_point = p;
		selected_oval.start(p);
	}

	public void extend_rotate(Point p) {
		selected_oval.rotate(p, start_point);
		repaint_all();
	}

	public void finish_rotate(Point p) {
		start_point = null;
	}

	public void clear() {
		JFrame frame = new JFrame();
		int ans = JOptionPane.showConfirmDialog(frame, "Are you sure you want to clear all?");
		if (ans == JOptionPane.YES_OPTION) {
			ovals.clear();
			repaint_all();
			operation_status = NONE;
		}
	}

	public int operation_status;
	public static final int NONE = 0;
	public static final int DRAWING = 1;
	public static final int SELECTING = 2;
	public static final int RESIZING = 3;
	public static final int MOVING = 4;
	public static final int ROTATING = 5;
	public Oval selected_oval = null;
	public int selected_handle = -1;

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
			if (operation_status == NONE) {
				for (Oval oval : ovals) {
					oval.end_selected();
				}
				selected_oval = null;
				operation_status = NONE;
				for (Oval oval : ovals) {
					if (oval.contains(p)) {
						selected_oval = oval;
						oval.start_selected();
						operation_status = SELECTING;
						repaint_all();
						break;
					}
				}
				repaint_all();
			}
			switch (operation_status) {
				case NONE:
					operation_status = DRAWING;
					start_stroke(p);
					break;
				case SELECTING:
					selected_handle = selected_oval.selected_handle(p);
					if (selected_handle >= 0) {
						operation_status = RESIZING;
						start_resize(p);
					} else if (selected_oval.contains(p)) {
						operation_status = MOVING;
						start_move(p);
					} else {
						operation_status = ROTATING;
						start_rotate(p);
					}
					break;
				default:
					break;
			}
		}

		public void mouseDragged(MouseEvent e) {
			Point p = e.getPoint();

			switch (operation_status) {
				case DRAWING:
					extend_stroke(p);
					break;
				case RESIZING:
					extend_resize(p);
					break;
				case MOVING:
					extend_move(p);
					break;
				case ROTATING:
					extend_rotate(p);
					break;
			}

		}

		public void mouseReleased(MouseEvent e) {
			Point p = e.getPoint();

			switch (operation_status) {
				case DRAWING:
					finish_stroke(p);
					operation_status = NONE;
					break;
				case RESIZING:
					finish_resize(p);
					operation_status = SELECTING;
					break;
				case MOVING:
					finish_move(p);
					operation_status = SELECTING;
					break;
				case ROTATING:
					finish_rotate(p);
					operation_status = SELECTING;
					break;
			}
		}

		public void mouseClicked(MouseEvent e) {
			Point p = e.getPoint();
			for (Oval oval : ovals) {
				oval.end_selected();
			}
			selected_oval = null;
			operation_status = NONE;
			for (Oval oval : ovals) {
				if (oval.contains(p)) {
					selected_oval = oval;
					oval.start_selected();
					operation_status = SELECTING;
					repaint_all();
					break;
				}
			}
			repaint_all();
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
		frame.setSize(1500, 1000);
		frame.setVisible(true);

	}
}
