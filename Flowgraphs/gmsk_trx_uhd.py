#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: GMSK, Ax.25, 9600, Streaming Mode
# Author: Zach Leffke, KJ4QLP
# Description: Half Duplex GMSK Ax.25 Transceiver, streaming processing
# GNU Radio version: v3.8.5.0-6-g57bd109d

from distutils.version import StrictVersion

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print("Warning: failed to XInitThreads()")

import os
import sys
sys.path.append(os.environ.get('GRC_HIER_PATH', os.path.expanduser('~/.grc_gnuradio')))

from PyQt5 import Qt
from gnuradio import qtgui
from gnuradio.filter import firdes
import sip
from ax25_gmsk_rx_hier import ax25_gmsk_rx_hier  # grc-generated hier_block
from ax25_gmsk_tx_hier import ax25_gmsk_tx_hier  # grc-generated hier_block
from datetime import datetime as dt; import string
from gnuradio import blocks
from gnuradio import filter
from gnuradio import gr
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import uhd
import time
from gnuradio.qtgui import Range, RangeWidget

from gnuradio import qtgui

class gmsk_trx_uhd(gr.top_block, Qt.QWidget):

    def __init__(self, radio_id='30CF9D2', rf_freq=401.12e6, rx_offset=250e3, sat_name='UPS1', tx_offset=250e3):
        gr.top_block.__init__(self, "GMSK, Ax.25, 9600, Streaming Mode")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("GMSK, Ax.25, 9600, Streaming Mode")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "gmsk_trx_uhd")

        try:
            if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                self.restoreGeometry(self.settings.value("geometry").toByteArray())
            else:
                self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        ##################################################
        # Parameters
        ##################################################
        self.radio_id = radio_id
        self.rf_freq = rf_freq
        self.rx_offset = rx_offset
        self.sat_name = sat_name
        self.tx_offset = tx_offset

        ##################################################
        # Variables
        ##################################################
        self.ts_str = ts_str = dt.strftime(dt.utcnow(), "%Y-%m-%dT%H:%M:%S.%fZ" )
        self.samp_rate = samp_rate = float(250000)
        self.fn = fn = "{:s}_{:s}_{:s}_{:s}k.fc32".format(sat_name, radio_id, ts_str, str(int(samp_rate/1e3)))
        self.tx_gain = tx_gain = 60
        self.tx_correct = tx_correct = 850
        self.rx_gain = rx_gain = 20
        self.rx_correct = rx_correct = 850
        self.interp = interp = 24
        self.fp = fp = "/captures/{:s}".format(fn)
        self.decim = decim = int(samp_rate/2000)
        self.bb_gain = bb_gain = .75

        ##################################################
        # Blocks
        ##################################################
        self._tx_gain_range = Range(0, 86, 1, 60, 200)
        self._tx_gain_win = RangeWidget(self._tx_gain_range, self.set_tx_gain, 'TX Gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._tx_gain_win, 4, 8, 1, 2)
        for r in range(4, 5):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(8, 10):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._tx_correct_range = Range(-10000, 10000, 1, 850, 200)
        self._tx_correct_win = RangeWidget(self._tx_correct_range, self.set_tx_correct, 'tx_correct', "counter_slider", float)
        self.top_grid_layout.addWidget(self._tx_correct_win, 5, 8, 1, 2)
        for r in range(5, 6):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(8, 10):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._rx_gain_range = Range(0, 86, 1, 20, 200)
        self._rx_gain_win = RangeWidget(self._rx_gain_range, self.set_rx_gain, 'RX Gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._rx_gain_win, 6, 0, 1, 2)
        for r in range(6, 7):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 2):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._rx_correct_range = Range(-10000, 10000, 1, 850, 200)
        self._rx_correct_win = RangeWidget(self._rx_correct_range, self.set_rx_correct, 'rx_correct', "counter_slider", float)
        self.top_grid_layout.addWidget(self._rx_correct_win, 6, 2, 1, 2)
        for r in range(6, 7):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(2, 4):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._bb_gain_range = Range(0, 1, .01, .75, 200)
        self._bb_gain_win = RangeWidget(self._bb_gain_range, self.set_bb_gain, 'TX bb_gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._bb_gain_win, 4, 10, 1, 2)
        for r in range(4, 5):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(10, 12):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("", "")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
        )
        self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(rf_freq+rx_correct, rx_offset), 0)
        self.uhd_usrp_source_0.set_gain(rx_gain, 0)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_now(uhd.time_spec(time.time()), uhd.ALL_MBOARDS)
        self.uhd_usrp_sink_0_0 = uhd.usrp_sink(
            ",".join(("", "")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
            '',
        )
        self.uhd_usrp_sink_0_0.set_center_freq(uhd.tune_request(rf_freq+tx_correct, tx_offset), 0)
        self.uhd_usrp_sink_0_0.set_gain(tx_gain, 0)
        self.uhd_usrp_sink_0_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0_0.set_samp_rate(250e3)
        self.uhd_usrp_sink_0_0.set_time_now(uhd.time_spec(time.time()), uhd.ALL_MBOARDS)
        self.rational_resampler_xxx_1 = filter.rational_resampler_ccc(
                interpolation=interp,
                decimation=decim,
                taps=None,
                fractional_bw=None)
        self.qtgui_waterfall_sink_x_0 = qtgui.waterfall_sink_c(
            2048, #size
            firdes.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            samp_rate / decim*interp, #bw
            "", #name
            1 #number of inputs
        )
        self.qtgui_waterfall_sink_x_0.set_update_time(0.010)
        self.qtgui_waterfall_sink_x_0.enable_grid(False)
        self.qtgui_waterfall_sink_x_0.enable_axis_labels(True)



        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        colors = [0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_waterfall_sink_x_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_waterfall_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_waterfall_sink_x_0.set_color_map(i, colors[i])
            self.qtgui_waterfall_sink_x_0.set_line_alpha(i, alphas[i])

        self.qtgui_waterfall_sink_x_0.set_intensity_range(-140, 10)

        self._qtgui_waterfall_sink_x_0_win = sip.wrapinstance(self.qtgui_waterfall_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_waterfall_sink_x_0_win, 4, 0, 2, 8)
        for r in range(4, 6):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 8):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.qtgui_freq_sink_x_1_0_0 = qtgui.freq_sink_c(
            2048, #size
            firdes.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            samp_rate/decim*interp, #bw
            "TX Spectrum", #name
            1
        )
        self.qtgui_freq_sink_x_1_0_0.set_update_time(0.010)
        self.qtgui_freq_sink_x_1_0_0.set_y_axis(-150, 0)
        self.qtgui_freq_sink_x_1_0_0.set_y_label('Relative Gain', 'dB')
        self.qtgui_freq_sink_x_1_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_1_0_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_1_0_0.enable_grid(True)
        self.qtgui_freq_sink_x_1_0_0.set_fft_average(1.0)
        self.qtgui_freq_sink_x_1_0_0.enable_axis_labels(True)
        self.qtgui_freq_sink_x_1_0_0.enable_control_panel(False)

        self.qtgui_freq_sink_x_1_0_0.disable_legend()


        labels = ['', '', '', '', '',
            '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
            "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(1):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_1_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_1_0_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_1_0_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_1_0_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_1_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_1_0_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_1_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_freq_sink_x_1_0_0_win, 0, 8, 4, 4)
        for r in range(0, 4):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(8, 12):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.qtgui_freq_sink_x_1_0 = qtgui.freq_sink_c(
            2048, #size
            firdes.WIN_BLACKMAN_hARRIS, #wintype
            0, #fc
            samp_rate / decim*interp, #bw
            "RX Spectrum", #name
            2
        )
        self.qtgui_freq_sink_x_1_0.set_update_time(0.010)
        self.qtgui_freq_sink_x_1_0.set_y_axis(-150, 0)
        self.qtgui_freq_sink_x_1_0.set_y_label('Relative Gain', 'dB')
        self.qtgui_freq_sink_x_1_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_1_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_1_0.enable_grid(True)
        self.qtgui_freq_sink_x_1_0.set_fft_average(1.0)
        self.qtgui_freq_sink_x_1_0.enable_axis_labels(True)
        self.qtgui_freq_sink_x_1_0.enable_control_panel(False)

        self.qtgui_freq_sink_x_1_0.disable_legend()


        labels = ['pre-d', 'agc_filt', '', '', '',
            '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
            1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
            "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0]

        for i in range(2):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_1_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_1_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_1_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_1_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_1_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_1_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_1_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_freq_sink_x_1_0_win, 0, 0, 4, 8)
        for r in range(0, 4):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 8):
            self.top_grid_layout.setColumnStretch(c, 1)
        self.blocks_socket_pdu_0_2 = blocks.socket_pdu('TCP_SERVER', '0.0.0.0', '8000', 1024, False)
        self.ax25_gmsk_tx_hier_0 = ax25_gmsk_tx_hier(
            bb_gain=0.75,
            bt=.5,
            samp_rate=250000,
        )
        self.ax25_gmsk_rx_hier_0 = ax25_gmsk_rx_hier(
            lpf_cutoff=12.5e3,
            lpf_trans=1e3,
            samp_rate=250000,
        )


        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.ax25_gmsk_rx_hier_0, 'kiss/ax25'), (self.blocks_socket_pdu_0_2, 'pdus'))
        self.msg_connect((self.blocks_socket_pdu_0_2, 'pdus'), (self.ax25_gmsk_tx_hier_0, 'kiss/ax25'))
        self.connect((self.ax25_gmsk_rx_hier_0, 1), (self.qtgui_freq_sink_x_1_0, 1))
        self.connect((self.ax25_gmsk_rx_hier_0, 0), (self.qtgui_freq_sink_x_1_0, 0))
        self.connect((self.ax25_gmsk_rx_hier_0, 1), (self.qtgui_waterfall_sink_x_0, 0))
        self.connect((self.ax25_gmsk_tx_hier_0, 0), (self.rational_resampler_xxx_1, 0))
        self.connect((self.ax25_gmsk_tx_hier_0, 0), (self.uhd_usrp_sink_0_0, 0))
        self.connect((self.rational_resampler_xxx_1, 0), (self.qtgui_freq_sink_x_1_0_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.ax25_gmsk_rx_hier_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "gmsk_trx_uhd")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_radio_id(self):
        return self.radio_id

    def set_radio_id(self, radio_id):
        self.radio_id = radio_id
        self.set_fn("{:s}_{:s}_{:s}_{:s}k.fc32".format(self.sat_name, self.radio_id, self.ts_str, str(int(self.samp_rate/1e3))))

    def get_rf_freq(self):
        return self.rf_freq

    def set_rf_freq(self, rf_freq):
        self.rf_freq = rf_freq
        self.uhd_usrp_sink_0_0.set_center_freq(uhd.tune_request(self.rf_freq+self.tx_correct, self.tx_offset), 0)
        self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(self.rf_freq+self.rx_correct, self.rx_offset), 0)

    def get_rx_offset(self):
        return self.rx_offset

    def set_rx_offset(self, rx_offset):
        self.rx_offset = rx_offset
        self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(self.rf_freq+self.rx_correct, self.rx_offset), 0)

    def get_sat_name(self):
        return self.sat_name

    def set_sat_name(self, sat_name):
        self.sat_name = sat_name
        self.set_fn("{:s}_{:s}_{:s}_{:s}k.fc32".format(self.sat_name, self.radio_id, self.ts_str, str(int(self.samp_rate/1e3))))

    def get_tx_offset(self):
        return self.tx_offset

    def set_tx_offset(self, tx_offset):
        self.tx_offset = tx_offset
        self.uhd_usrp_sink_0_0.set_center_freq(uhd.tune_request(self.rf_freq+self.tx_correct, self.tx_offset), 0)

    def get_ts_str(self):
        return self.ts_str

    def set_ts_str(self, ts_str):
        self.ts_str = ts_str
        self.set_fn("{:s}_{:s}_{:s}_{:s}k.fc32".format(self.sat_name, self.radio_id, self.ts_str, str(int(self.samp_rate/1e3))))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_decim(int(self.samp_rate/2000))
        self.set_fn("{:s}_{:s}_{:s}_{:s}k.fc32".format(self.sat_name, self.radio_id, self.ts_str, str(int(self.samp_rate/1e3))))
        self.qtgui_freq_sink_x_1_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)
        self.qtgui_freq_sink_x_1_0_0.set_frequency_range(0, self.samp_rate/self.decim*self.interp)
        self.qtgui_waterfall_sink_x_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)

    def get_fn(self):
        return self.fn

    def set_fn(self, fn):
        self.fn = fn
        self.set_fp("/captures/{:s}".format(self.fn))

    def get_tx_gain(self):
        return self.tx_gain

    def set_tx_gain(self, tx_gain):
        self.tx_gain = tx_gain
        self.uhd_usrp_sink_0_0.set_gain(self.tx_gain, 0)

    def get_tx_correct(self):
        return self.tx_correct

    def set_tx_correct(self, tx_correct):
        self.tx_correct = tx_correct
        self.uhd_usrp_sink_0_0.set_center_freq(uhd.tune_request(self.rf_freq+self.tx_correct, self.tx_offset), 0)

    def get_rx_gain(self):
        return self.rx_gain

    def set_rx_gain(self, rx_gain):
        self.rx_gain = rx_gain
        self.uhd_usrp_source_0.set_gain(self.rx_gain, 0)

    def get_rx_correct(self):
        return self.rx_correct

    def set_rx_correct(self, rx_correct):
        self.rx_correct = rx_correct
        self.uhd_usrp_source_0.set_center_freq(uhd.tune_request(self.rf_freq+self.rx_correct, self.rx_offset), 0)

    def get_interp(self):
        return self.interp

    def set_interp(self, interp):
        self.interp = interp
        self.qtgui_freq_sink_x_1_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)
        self.qtgui_freq_sink_x_1_0_0.set_frequency_range(0, self.samp_rate/self.decim*self.interp)
        self.qtgui_waterfall_sink_x_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)

    def get_fp(self):
        return self.fp

    def set_fp(self, fp):
        self.fp = fp

    def get_decim(self):
        return self.decim

    def set_decim(self, decim):
        self.decim = decim
        self.qtgui_freq_sink_x_1_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)
        self.qtgui_freq_sink_x_1_0_0.set_frequency_range(0, self.samp_rate/self.decim*self.interp)
        self.qtgui_waterfall_sink_x_0.set_frequency_range(0, self.samp_rate / self.decim*self.interp)

    def get_bb_gain(self):
        return self.bb_gain

    def set_bb_gain(self, bb_gain):
        self.bb_gain = bb_gain




def argument_parser():
    description = 'Half Duplex GMSK Ax.25 Transceiver, streaming processing'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        "--radio-id", dest="radio_id", type=str, default='30CF9D2',
        help="Set radio_id [default=%(default)r]")
    parser.add_argument(
        "--rf-freq", dest="rf_freq", type=eng_float, default="401.12M",
        help="Set rf_freq [default=%(default)r]")
    parser.add_argument(
        "--rx-offset", dest="rx_offset", type=eng_float, default="250.0k",
        help="Set rx_offset [default=%(default)r]")
    parser.add_argument(
        "--sat-name", dest="sat_name", type=str, default='UPS1',
        help="Set sat_name [default=%(default)r]")
    parser.add_argument(
        "--tx-offset", dest="tx_offset", type=eng_float, default="250.0k",
        help="Set tx_offset [default=%(default)r]")
    return parser


def main(top_block_cls=gmsk_trx_uhd, options=None):
    if options is None:
        options = argument_parser().parse_args()

    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls(radio_id=options.radio_id, rf_freq=options.rf_freq, rx_offset=options.rx_offset, sat_name=options.sat_name, tx_offset=options.tx_offset)

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    def quitting():
        tb.stop()
        tb.wait()

    qapp.aboutToQuit.connect(quitting)
    qapp.exec_()

if __name__ == '__main__':
    main()
