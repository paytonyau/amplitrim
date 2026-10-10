# Contributing

Bug reports, questions and improvements are welcome. Use the repository's [issues](https://github.com/paytonyau/amplitrim/issues) for all of them.

## Reporting a problem or asking a question

Please include:

- the output of `amplitrim --version`, your operating system, and the version of vsearch (`vsearch --version`);
- the command you ran, and what you expected;
- `REPORT.md` and `summary.json` from the run, if it got that far. They contain no reads. Please do not attach your FASTQ files.

If the tool recommends something that looks wrong for your data, that is worth reporting even if it did not fail.

## Changing the code

1. Fork the repository and install it with the test extras: `pip install -e ".[test]"`. vsearch must be on your PATH for the tests that run it; the others run without it.
2. Run `python -m pytest tests -q` before and after your change.
3. Add a test for what you change, and keep `README.md` and `docs/usage.md` in step with the behaviour.
4. Open a pull request that says what it changes and why.

People use AmpliTrim's recommendations to choose settings. A change to **what it recommends** (not to how it is run or reported) therefore needs evidence on data, in the way the [validation page](https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md) gives it. Write the criteria down before the run, and report the outcome whether or not it holds. The validation scripts are not in this repository; they are in the [validation archive](https://github.com/paytonyau/amplitrim-validation).

## Conduct

Be courteous and assume good faith. Criticism of the code and of the validation is welcome and is what the tests and the written criteria are for.
