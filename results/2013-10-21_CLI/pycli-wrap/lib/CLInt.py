#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import msvcrt
import time

import Common

class CLImain(Common.CLIcommon):

    def Prompt(self):
        prompt = "\n" + self.PromptStr
        sys.stdout.write(prompt)
        sys.stdout.flush()
    
    def run(self):
        input_text = ""
        sw = 0
        self.InitialTasks()
        self.Prompt()
        try:
            while True:
            
                if msvcrt.kbhit():
                    readit = msvcrt.getch()
                    sys.stdout.write(readit)
                    input_text += readit
                    
                    if ord(readit) in self.EnterKeys():
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
