package board;

//	base code for simple applications
//  Copyright (c) 2002 Takeo Igarashi 

//  basic button class 
//	use with ButtonsInterface.java

import java.awt.event.*;
import java.awt.*;
import java.util.*;

interface ButtonsInterface
{
	public void buttonPressed(String label);
	public void initializeButtons(Buttons buttons);
}

public class Buttons extends Panel {


    ButtonsInterface target;
    TextField filename;

	static public Font font = new Font("Helvetica", Font.PLAIN, 12);

	
    public Buttons(ButtonsInterface target) {
		this.target = target;
		setLayout(new FlowLayout());
		setBackground(Color.lightGray);

	
		target.initializeButtons(this);
		
    }
	public void buttonPressed(String name){
		target.buttonPressed(name);
	}

	
	public void addButton(String label){
		add(new MyButton(label));
	}
	
	
	
	
	
	
	
	
    public void paint(Graphics g) {
		Rectangle r = getBounds();

		g.setColor(Color.lightGray);
		g.draw3DRect(0, 0, r.width, r.height, false);
    }

	
	class MyButton extends Button {
		MyButton(String label){
			setFont(Buttons.font);
			setLabel(label);
			addMouseListener(
				new MouseAdapter() {
      				public void mousePressed(MouseEvent e) {
						buttonPressed(getLabel());
					}
				}
			);
		}
	}

}


