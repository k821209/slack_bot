import os
import time
from slackclient import SlackClient
import kang
import pandas as pd 
import primer3



dicFA    = kang.Fasta2dic('../References/Athaliana/Araport11_genes.201606.cdna.fasta')
dicFA_func = kang.Fasta2dic_all('../References/Athaliana/Araport11_genes.201606.cdna.fasta')
genelist = dicFA.keys()
keys   = [x.split('|')[0].replace('>','').strip() for x in dicFA_func.keys()]
values = [x.split('|')[1].strip() for x in dicFA_func.keys()]
dicG2F = dict(zip(keys,values))


# starterbot's ID as an environment variable
BOT_ID = os.environ.get("BOT_ID")

# constants
AT_BOT = "<@" + BOT_ID + ">"
COMMAND_list = ["seq","func","primer"]

# instantiate Slack & Twilio clients
slack_client = SlackClient(os.environ.get('SLACK_BOT_TOKEN'))

def get_opt_cloningprimer_pair(seq):
    oligo_calc = primer3.thermoanalysis.ThermoAnalysis(mv_conc=20, dv_conc=1.5, dntp_conc=0.8, dna_conc=50, max_nn_length=60)
    # Change if you have personal condition for PCR 
    ## mv_conc   : The millimolar (mM) concentration of monovalent salt cations (usually KCl) in the PCR.
    ## dv_conc   : The millimolar concentration of divalent salt cations (usually MgCl^(2+)) in the PCR.
    ## dntp_conc : The millimolar concentration of the sum of all deoxyribonucleotide triphosphates.
    ## dna_conc  : A value to use as nanomolar (nM) concentration of each annealing oligo over the course the PCR. 
    ## max_nn_length : longest seq length for primer Tm analysis
    primerF_list = []
    primerR_list = []
    for i in range(16,25):
        primerF = seq[0:i]
        #primerR = kang.rev_comp(seq[-i-3:-3]) # stop codon excluded
        primerR = kang.rev_comp(seq[-i:]) # stop codon included
        if 45 <= oligo_calc.calcTm(primerF) <= 63:
            primerF_list.append(primerF)
        else:
            pass
            #print (primerF, oligo_calc.calcTm(primerF) )
        if 45 <= oligo_calc.calcTm(primerR) <= 63:
            primerR_list.append(primerR)

    primerF_list_Tm = [oligo_calc.calcTm(x) for x in primerF_list]
    primerR_list_Tm = [oligo_calc.calcTm(x) for x in primerR_list]

    return primerF_list,primerR_list,primerF_list_Tm,primerR_list_Tm



def handle_command(command, channel):
    """
        Receives commands directed at the bot and determines if they
        are valid commands. If so, then acts on the commands. If not,
        returns back what it needs for clarification.
    """
    response = \
    '''
	Use commands 
        "seq %s"
        "func %s"
        "primer %s"
    '''%(genelist[0],genelist[1],genelist[2])
    if command.startswith(COMMAND_list[0]): # seq
        seq = kang.txtwrap(dicFA[command.split()[1].strip().upper()],50)
        response = '\n'.join(seq)
    elif command.startswith(COMMAND_list[1]): # func
        response = dicG2F[command.split()[1].strip().upper()]
    elif command.startswith(COMMAND_list[2]): # primer
        seq = dicFA[command.split()[1].strip().upper()]
        primerF_list,primerR_list,primerF_list_Tm,primerR_list_Tm = get_opt_cloningprimer_pair(seq)
        primerF = ["%s (%f)"%(x,y) for x,y in zip(primerF_list,primerF_list_Tm)]
        primerR = ["%s (%f)"%(x,y) for x,y in zip(primerR_list,primerR_list_Tm)]
        print primerF,primerR
        response = """-Forward primers-
%s
-Reverse primers-
%s
"""%('\n'.join(primerF),'\n'.join(primerR))
    slack_client.api_call("chat.postMessage", channel=channel,
                          text=response, as_user=True)


def parse_slack_output(slack_rtm_output):
    """
        The Slack Real Time Messaging API is an events firehose.
        this parsing function returns None unless a message is
        directed at the Bot, based on its ID.
    """
    output_list = slack_rtm_output
    if output_list and len(output_list) > 0:
        for output in output_list:
            if output and 'text' in output and AT_BOT in output['text']:
                # return text after the @ mention, whitespace removed
                return output['text'].split(AT_BOT)[1].strip().lower(), \
                       output['channel']
    return None, None

if __name__ == "__main__":
    READ_WEBSOCKET_DELAY = 1 # 1 second delay between reading from firehose
    if slack_client.rtm_connect():
        print("StarterBot connected and running!")
        while True:
            command, channel = parse_slack_output(slack_client.rtm_read())
            if command and channel:
                handle_command(command, channel)
            time.sleep(READ_WEBSOCKET_DELAY)
    else:
        print("Connection failed. Invalid Slack token or bot ID?")


