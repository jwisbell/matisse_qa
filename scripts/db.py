import sqlite3


# Create and connect to the SQLite database
def create_database(db_name):
    conn = sqlite3.connect(f"{db_name}")
    cursor = conn.cursor()

    # Create the table with a UNIQUE constraint on OBS_TIME
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Observations (
            TARGET TEXT,
            OBS_TIME TEXT UNIQUE,
            TAU REAL
        )
    """)
    conn.commit()
    conn.close()


# Function to insert a new observation only if OBS_TIME is not in the database
def insert_observation(target, obs_time, tau0, db):
    conn = sqlite3.connect(f"{db}")
    cursor = conn.cursor()

    # Attempt to insert the observation
    try:
        cursor.execute(
            """
            INSERT INTO Observations (TARGET, OBS_TIME, TAU)
            VALUES (?, ?,?)
        """,
            (target, obs_time, round(tau0, 2)),
        )
        conn.commit()
        print("Observation added successfully.")
    except sqlite3.IntegrityError:
        print("Observation with this OBS_TIME already exists. No new entry added.")

    conn.close()


def get_obs(db_name):
    with sqlite3.connect(db_name) as conn:
        cur = conn.cursor()
        cur.execute("SELECT * FROM Observations")
        rows = list(cur.fetchall())
        for row in rows:
            print(row)


if __name__ == "__main__":
    from sys import argv

    script, db_name = argv
    get_obs(db_name)
