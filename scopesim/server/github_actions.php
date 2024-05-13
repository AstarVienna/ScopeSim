<?php
  $owner = $_GET["owner"] ?? "AstarVienna";
  $repo = $_GET["repo"] ?? "ScopeSim";
  $workflow_file = $_GET["workflow_file"] ?? "tests.yml";
  $branch = $_GET["branch"] ?? "main";

  $git_run_url = "https://api.github.com/repos/AstarVienna/".$owner."/".$repo."/workflows/".$workflow_file."/runs?branch=".$branch."&per_page=1";
  $data = json_decode(file_get_contents($git_run_url), true);

  $python_version = "3.7";
  $run_status = "passing";
  $colour = "brightgreen";

  echo $data["workflow_runs"][0]["id"]

  # $git_job_url = "https://api.github.com/repos/AstarVienna/ScopeSim/actions/runs/{run_id}/jobs"
  # $data = json_decode(file_get_contents($git_run_url), true);

  # $badge_url = "https://img.shields.io/badge/Python_" . $python_version . "-" . $run_status . "-" . $colour;
  # $badge = file_get_contents($git_run_url), true);

  # echo header('Location: ' . $badge_url);
  
?>