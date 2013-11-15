#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, select, tty, termios

import Common

#Thanks to Graham King --> http://www.darkcoding.net/software/non-blocking-console-io-is-not-possible/
class CLImain(Common.CLIcommon):

    def isData(self):
	return select.select([sys.stdin], [], [], 0) == ([sys.stdin], [], [])
    
    def Prompt(self):
	sys.stdout.write(self.PromptStr)
	sys.stdout.flush()
    
    def run(self):
	old_settings = termios.tcgetattr(sys.stdin)
	input_text = ""
	sw = 0
        self.InitialTasks()
        self.Prompt()
	
	try:
	    tty.setcbreak(sys.stdin.fileno())
	    while 1:
		if self.isData():
		    
		    readit = sys.stdin.read(1)
		    sys.stdout.write(readit)
		    input_text += readit
		    
		    if readit == '\x1b':         # x1b is ESC
			break
		    
		    elif ord(readit) in self.EnterKeys():
			self.PressEnterTasks()
			if input_text.strip() in self.ExitWords() or ord(readit) in self.ExitKeys():
			    self.CleanUpTasks()

			if input_text.strip() != "":
			    self.ProcessCLI(input_text)
			    
			    
			self.Prompt()
			input_text = ""
		
		self.InfiniteLoopTasks()
	    
	except KeyboardInterrupt:
	    self.CleanUpTasks()
		    
	finally:
		termios.tcsetattr(sys.stdin, termios.TCSADRAIN, old_settings)
