import jinja2
import os


def render_jinja_by_params(template_file, rendered_file, params):
    template_dir = os.path.dirname(template_file)
    template_name = os.path.basename(template_file)

    template_loader = jinja2.FileSystemLoader(searchpath=template_dir)
    template_env = jinja2.Environment(loader=template_loader)
    template_env.globals.update(zip=zip)
    template = template_env.get_template(template_name)
    rendered_script = template.render(
        **params
    )

    with open(rendered_file, 'w') as out_file:
        out_file.write(rendered_script)


if __name__ == "__main__":
    render_jinja_by_params(
        str(snakemake.input[0]),
        str(snakemake.output[0]),
        snakemake.params
    )
