import os
from datetime import datetime

__version__ = '1.3.0'
__version_full__ = __version__

def get_path():
    """
    Shortcut for users whose theme is next to their conf.py.
    """
    # Theme directory is defined as our parent directory
    return os.path.abspath(os.path.dirname(os.path.dirname(__file__)))


def update_context(app, pagename, templatename, context, doctree):
    context['nltk_theme_version'] = __version_full__

def setup(app):
    # add_html_theme is new in Sphinx 1.6+
    if hasattr(app, 'add_html_theme'):
        theme_path = os.path.abspath(os.path.dirname(__file__))
        app.add_html_theme('nltk_theme', theme_path)
    app.connect('html-page-context', update_context)
    return {'version': '0.3.0',
            'parallel_read_safe': True}
