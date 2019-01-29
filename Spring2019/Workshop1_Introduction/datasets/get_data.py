#!/usr/bin/env python
from selenium.webdriver import Firefox
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support.ui import Select
opts = Options()
#opts.set_headless()
#assert opts.headless
browser = Firefox(options=opts)

browser.get("http://hbcponline.com/faces/contents/InteractiveTable.xhtml?search=false&tableId=6")
export_data_button = browser.find_element_by_id('leftForm:sm_export')
export_data_button.click()
#rowselect_radio = browser.find_element_by_id('tblExportForm:rowsel:0')

#filetype_selector = browser.find_element_by_id('tblExportForm:fileTypeMenu_input')
#filetype_selector.select_by_value('csv')

select = Select(browser.find_element_by_id("tblExportForm:fileTypeMenu_input"))
print(select.options)
browser.close()
                



#export_form = browser.find_element_by_id('tblExportDialog')

#browser.close()
