import os
import socket
import click
from seqnado_web.app import app

def get_ip():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't even have to be reachable
        s.connect(('10.255.255.255', 1))
        IP = s.getsockname()[0]
    except Exception:
        IP = '127.0.0.1'
    finally:
        s.close()
    return IP

@click.command()
@click.option("--port", default=5001, help="Port to run the server on")
@click.option("--host", default="0.0.0.0", help="Host to run the server on")
@click.option("--debug/--no-debug", default=False, help="Run in debug mode")
def main(port, host, debug):
    """Run the SeqNado Web Interface."""
    hostname = socket.gethostname()
    user = os.getenv("USER", "user")
    
    print("\n" + "="*60)
    print(f"🚀 SeqNado Web App starting on port {port}...")
    print("="*60)
    print("\nTo access this interface from your local machine, run:")
    print(f"\n    ssh -L {port}:localhost:{port} {user}@{hostname}\n")
    print(f"Then open your browser to: http://localhost:{port}")
    print("="*60 + "\n")
    
    if debug:
        print("Starting Flask development server...")
        app.run(host=host, port=port, debug=True)
    else:
        import shutil
        import subprocess
        
        if not shutil.which("gunicorn"):
            print("WARNING: 'gunicorn' not found. Falling back to Flask development server.")
            print("To remove this warning, install gunicorn: pip install gunicorn")
            app.run(host=host, port=port, debug=False)
        else:
            print("Starting Gunicorn production server...")
            # Run gunicorn
            cmd = [
                "gunicorn",
                "-w", "4",
                "-b", f"{host}:{port}",
                "--access-logfile", "-", 
                "--error-logfile", "-",
                "seqnado_web.app:app"
            ]
            try:
                subprocess.run(cmd, check=True)
            except KeyboardInterrupt:
                pass
            except Exception as e:
                print(f"Error running gunicorn: {e}")
                app.run(host=host, port=port, debug=False)

if __name__ == "__main__":
    main()
