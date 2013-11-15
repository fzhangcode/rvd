#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import time


class CLIcommon(object):
    
    def __init__(self):
        self.PromptStr = "CLI> "
        self.InitialTasks = self.func_InitialTasks
        self.EnterKeys = self.func_EnterKeys
        self.PressEnterTasks = self.func_PressEnterTasks
        self.ExitWords = self.func_ExitWords
        self.ExitKeys = self.func_ExitKeys
        self.CleanUpTasks = self.func_CleanUpTasks
        self.InfiniteLoopTasks = self.func_InfiniteLoopTasks
        self.ProcessCLI_dict = {
            'help': self.help,
            'dummy': self.dummy,
        }
        self.ProcessCLI = self.func_ProcessCLI
        

    def configure(self, **kargs):
        for name,value in kargs.items():
            isExist = hasattr(self, name)
            if isExist:
                setattr(self, name, value)
            else:
                print "ABORT EXECUTION: No existe funcion ", name
                sys.exit(1)
              
    def func_InitialTasks(self):
        self.func_Banner()
        
    def func_Banner(self):
        print """

###################################
# Welcome to CLI Python Prototype #
###################################
    
        """

    def func_EnterKeys(self, enter_keys = [10, 13, 4, 3]):
        return enter_keys
    
    def func_PressEnterTasks(self):
        mytime = time.localtime()
        mystrtime = time.strftime("%d-%m-%Y %H:%M", mytime)
        print "Hora actual: " + mystrtime
        
    def func_ExitWords(self, exit_words = ['exit', 'salir', 'quit']):
        return exit_words
    
    def func_ExitKeys(self, exit_keys = [4, 3]):
        return exit_keys

    def func_CleanUpTasks(self):
        print "Exiting ..."
        sys.exit(0)

    def func_InfiniteLoopTasks(self):
        pass

    def func_ProcessCLI(self, command):
        """
        """

        command = command.strip()
        if command in self.ProcessCLI_dict:
            self.ProcessCLI_dict[command]()
            
        elif command == "this.valid_cmds":
            return self.ProcessCLI_dict
        
        else:
            print "Error: Comando incorrecto. "

    def help(self):
        """Muestra la ayuda sobre los comandos disponibles."""
    
        helptext = """
        
    ### Comandos disponibles:\n"""
        valid_cmds = self.ProcessCLI("this.valid_cmds")
    
        for key, value in valid_cmds.items():
            helptext += "\n\t"
            helptext += key + " - " + value.__doc__
           
        helptext += ""    
        print helptext
    
    def dummy(self):
        """Comando para pruebas."""
    
        print """
            Dummy rulez!
        """
        
