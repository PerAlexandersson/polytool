//! Generate a curated family and prepare a reproducible OEIS draft.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let catalog = polytool::families::catalog()?;
    let family = catalog
        .iter()
        .find(|f| f.id == "PT-F0001")
        .ok_or("missing family")?;
    for row in family.generate_rows(6)? {
        println!("{}", row.join(", "));
    }
    let draft = polytool::families::export_oeis(family, 12)?;
    println!("{}", draft.blocks);
    Ok(())
}
