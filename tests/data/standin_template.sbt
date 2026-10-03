Submit-block ::= {
  contact {
    contact {
      name name {
        last "Submitter",
        first "Example",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "Example University",
        div "Example Department",
        city "Example City",
        sub "",
        country "Vietnam",
        street "1 Example Street",
        email "example@example.org",
        postal-code "000000"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Submitter",
            first "Example",
            middle "",
            initials "",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "Example University",
        div "Example Department",
        city "Example City",
        sub "",
        country "Vietnam",
        street "1 Example Street",
        postal-code "000000"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Submitter",
              first "Example",
              middle "",
              initials "",
              suffix "",
              title ""
            }
          }
        }
      },
      title "STAND-IN submission template - NOT FOR SUBMISSION"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "This template is a test fixture: the names and the address are placeholders."
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
