from selenium import webdriver
from selenium.webdriver.common.keys import Keys

driver = webdriver.Firefox()
driver.get("https://services.healthtech.dtu.dk/service.php?TMHMM-2.0") #URL of wanted page
'''
assert "Python" in driver.title
elem = driver.find_element_by_name("q")
elem.clear()
elem.send_keys("pycon")
elem.send_keys(Keys.RETURN)
assert "No results found." not in driver.page_source
driver.close()
'''

<inpout type = "text" name = "seq" id = "seq_id"/> 




