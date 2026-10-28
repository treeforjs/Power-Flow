"""Verify the actual compressed source-acquisition packet without source requests."""
import json
from packet_common import WORK, contract, verify_archive, write_new


def main():
    result = verify_archive(contract())
    write_new(WORK/'delivery_verification.json', result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
