"""
This file contains most of the configuration variables that your app needs.
"""

DEBUG = False  # Turns on debugging features in Flask
MAIL_FROM_EMAIL = "raphamontana@gmail.com"  # For use in application emails
BCRYPT_LOG_ROUNDS = 12  # Configuration for the Flask-Bcrypt extension

# If you're using Flask-Bcrypt to hash user passwords, you'll need to specify the number of "rounds" that the algorithm executes in hashing a password.
# If you aren't using Flask-Bcrypt, you should probably start.
# The more rounds used to hash a password, the longer it'll take for an attacker to guess a password given the hash.
# The number of rounds should increase over time as computing power increases.

