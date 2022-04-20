
# from lamp.start_lamp import LAMP
# from design.search_primers import Design


# class LoadingThread(QThread):

#     progress_update = pyqtSignal(int)

#     def __init__(self):
#         QThread.__init__(self)
#         self.current_class = None
#         self.loading_progress = 0


#     def __del__(self):
#         self.wait()
    

#     def run(self) -> None:
#         while True:
#             self.progress_update.emit()
#             time.sleep(1)