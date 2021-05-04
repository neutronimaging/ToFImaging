from JupyterNotebooks.code.parent import Parent


class GuiHandler(Parent):

    def moving_average_size_radioButton_clicked(self):
        state_default_kernel = self.parent.ui.kernel_size_default_radioButton.isChecked()
        self.parent.ui.kernel_size_default_label.setEnabled(state_default_kernel)
        list_ui = [self.parent.ui.kernel_size_custom_y_label,
                   self.parent.ui.kernel_size_custom_y_lineEdit,
                   self.parent.ui.kernel_size_custom_x_label,
                   self.parent.ui.kernel_size_custom_x_lineEdit,
                   self.parent.ui.kernel_size_custom_lambda_label,
                   self.parent.ui.kernel_size_custom_lambda_lineEdit]
        for _ui in list_ui:
            _ui.setEnabled(not state_default_kernel)
