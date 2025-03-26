#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/3/21 12:28
# @Author  : miko
# @File    : crawler.py
# @Software: PyCharm
# @Site:

import os
import time
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By

from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.action_chains import ActionChains

# Set Downloader
chrome_options = webdriver.ChromeOptions()
download_dir = os.path.join(os.getcwd(), "../downloads")
os.makedirs(download_dir, exist_ok=True)

prefs = {
    "download.default_directory": download_dir,
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True
}
chrome_options.add_argument("--window-size=1920,1080")
# chrome_options.add_argument("--headless")
chrome_options.add_experimental_option("prefs", prefs)

# Example SMILES: COCCN1C=C(C=N1)C2=CC=CC3=C2C=CC=C3N C1=CC=C(C=C1)N2C=C(C3=CC=CC=C32)/C=C(\CN)/C(=O)O

# Init Driver
driver = webdriver.Chrome(options=chrome_options)
driver.get("http://127.0.0.1:7800/")

time.sleep(5)

actions = ActionChains(driver)
actions.move_by_offset(200, 300).click().perform()

time.sleep(1)

actions = ActionChains(driver)
actions.key_down(Keys.CONTROL).send_keys('v').key_up(Keys.CONTROL).perform()

time.sleep(1)

textarea_element = driver.find_element(By.TAG_NAME, "textarea")
textarea_element.send_keys('COCCN1C=C(C=N1)C2=CC=CC3=C2C=CC=C3N')
textarea_element.send_keys(Keys.RETURN)

input()

accept_button = driver.find_element(By.XPATH, "//div[contains(@class,'gwt-DecoratedPopupPanel')]//button[text()='Accept']")
print(accept_button)
accept_button.click()

# actions = ActionChains(driver)
# actions.move_by_offset(1341.23, 665).click().perform()

input()