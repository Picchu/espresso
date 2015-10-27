#from __future__ import absolute_import

from celery import Celery

# configuration for celery app using rabbitmq as both
# broker and backend (necessary to retrieve results from
# tasks)
my_app = Celery('lib',
             broker='amqp://',
             backend='amqp://',
             include=['lib.tasks_app'])

# Optional configuration, see the application user guide.
my_app.conf.update(
    CELERY_TASK_RESULT_EXPIRES=3600,
)

if __name__ == '__main__':
    my_app.start()


# in order to start the worker, run this command in
# another terminal in the folder above this one
#celery worker --app=lib.celery_app:my_app -l info

# TO-DO: call the worker automatically within python program
