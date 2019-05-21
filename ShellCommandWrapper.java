/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package mpjtest;
import java.io.IOException;
import java.util.*;


/**
 *
 * @author ahsan
 */
public class ShellCommandWrapper
{
  public static void shellCommandExecutor(String[] args)
  throws IOException, InterruptedException
  {
    // you need a shell to execute a command pipeline
    List<String> commands = new ArrayList<String>();
    for (int i=0; i< args.length; i++){
        commands.add(args[i]);
        System.out.println(args[i]);
    }
    SystemCommandExecutor commandExecutor = new SystemCommandExecutor(commands);
    int result = commandExecutor.executeCommand();
 
    StringBuilder stdout = commandExecutor.getStandardOutputFromCommand();
    StringBuilder stderr = commandExecutor.getStandardErrorFromCommand();
 
    System.out.println("STDOUT");
    System.out.println(stdout);
 
    System.out.println("STDERR");
    System.out.println(stderr);
  }
}
