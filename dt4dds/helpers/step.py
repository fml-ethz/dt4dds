import time
import gc
import pprint

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class Step():


    def __init__(self):
        
        self._pre_process_hooks = []
        self._post_process_hooks = []
        
        self._pre_loop_hooks = []
        self._post_loop_hooks = []

        self._start_time = None


    def __repr__(self):
        return f"{type(self).__name__}()"


    def process(self, *args, **kwargs):
        raise NotImplementedError(f"Method process has to be implemented in the subclass of Step.")



    def add_post_process_hook(self, callable):
        """"""
        self._post_process_hooks.append(callable)


    def add_pre_process_hook(self, callable):
        """"""
        self._pre_process_hooks.append(callable)

    
    def _run_post_process_hooks(self, *args, **kwargs):
        for callable in self._post_process_hooks:
            try:
                callable(*args, **kwargs)
                logger.info(f"Successfully executed post process hook '{callable.__name__}' for '{self}'.")
            except Exception:
                logger.exception(f"Unable to process post process hook '{callable.__name__}' for '{self}'.")

        logger.info(f"Step {self} took {time.time()-self._start_time:.1f}s to complete.")
        gc.collect()


    def _run_pre_process_hooks(self, *args, **kwargs):
        gc.collect()
        self._start_time = time.time()
        logger.info(f"Running {self}.")

        # print settings info if available
        if hasattr(self, 'settings') and logger.level <= logging.INFO:
            settings_dict = pprint.pformat(self.settings, compact=True, width=10000000)
            logger.debug(f"Settings: {settings_dict}")

        for callable in self._pre_process_hooks:
            try:
                callable(*args, **kwargs)
                logger.info(f"Successfully executed pre process hook '{callable.__name__}' for '{self}'.")
            except Exception:
                logger.exception(f"Unable to process pre process hook '{callable.__name__}' for '{self}'.")



    def add_post_loop_hook(self, callable):
        """"""
        self._post_loop_hooks.append(callable)


    def add_pre_loop_call(self, callable):
        """"""
        self._pre_loop_hooks.append(callable)

    
    def _run_post_loop_hooks(self, *args, **kwargs):
        for callable in self._post_loop_hooks:
            try:
                callable(*args, **kwargs)
                logger.info(f"Successfully executed post loop hook '{callable.__name__}' for '{self}'.")
            except Exception:
                logger.exception(f"Unable to process post loop hook '{callable.__name__}' for '{self}'.")


    def _run_pre_loop_hooks(self, *args, **kwargs):
        for callable in self._pre_loop_hooks:
            try:
                callable(*args, **kwargs)
                logger.info(f"Successfully executed pre loop hook '{callable.__name__}' for '{self}'.")
            except Exception:
                logger.exception(f"Unable to process pre loop hook '{callable.__name__}' for '{self}'.")

