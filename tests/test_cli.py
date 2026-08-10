"""Tests for the hydrolib CLI."""

from click.testing import CliRunner

from hydrolib.cli import cli


class TestSkewCommand:
    def test_known_state_prints_value_and_citation(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["skew", "Vermont"])
        assert result.exit_code == 0
        assert "0.440" in result.output
        assert "Olson, S.A., 2014" in result.output

    def test_accepts_abbreviation(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["skew", "VT"])
        assert result.exit_code == 0
        assert "0.440" in result.output

    def test_unlisted_state_prints_nationwide_fallback(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["skew", "California"])
        assert result.exit_code == 0
        assert "-0.302" in result.output
        assert "Bulletin 17B" in result.output

    def test_invalid_state_exits_with_error(self):
        runner = CliRunner()
        result = runner.invoke(cli, ["skew", "Not A State"])
        assert result.exit_code != 0
        assert "Unrecognized state" in result.output
