import os
import sys
from google import genai

# Pass the gemini api key
api_key = os.environ.get('GEMINI_API_KEY')
client = genai.Client(api_key=api_key)

prompt = sys.stdin.read() # Expects a git diff

response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents="Create a concise summary in markdown format of the following git diff of a pull request: \n" + prompt
)

summary = response.text

# Strip the response from a surrounding markdown code block
if summary.startswith("```markdown"):
    summary = summary.removeprefix("```markdown").strip()
elif summary.startswith("```"):
    summary = summary.removeprefix("```").strip()

if summary.endswith("```"):
    summary = summary.removesuffix("```").strip()

print(summary)