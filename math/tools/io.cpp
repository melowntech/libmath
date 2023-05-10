#include "service/cmdline.hpp"
#include "utility/buildsys.hpp"
#include "utility/enum-io.hpp"

#include "math/geometry_core.hpp"

namespace po = boost::program_options;

namespace {

UTILITY_GENERATE_ENUM(Type,
                      ((Size2))
                      ((Size2i))
                      ((Size2ll))
                      ((Size2r))
                      ((Size2f))

                      ((Size3))
                      ((Size3i))
                      ((Size3ll))
                      ((Size3r))
                      ((Size3f))

                      ((Extents2))
                      ((Extents2i))
                      ((Extents2ll))
                      ((Extents2f))

                      ((Extents3))
                      ((Extents3i))
                      ((Extents3ll))
                      ((Extents3f))

                      ((Viewport2))
                      ((Viewport2i))
                      ((Viewport2f))
                      )

struct LexicalCast {
    Type type;
    std::string value;

    using optional = std::optional<LexicalCast>;

    LexicalCast(Type type, std::string value)
        : type(type), value(std::move(value))
    {}
};

class MathIo : public service::Cmdline
{
public:
    MathIo()
        : Cmdline("math.io", BUILD_TARGET_VERSION
                  , service::DISABLE_EXCESSIVE_LOGGING)
    {
    }

private:
    virtual void configuration(po::options_description &cmdline
                               , po::options_description &config
                               , po::positional_options_description &pd)
        override;

    virtual void configure(const po::variables_map &vars)
        override;

    virtual bool help(std::ostream &out, const std::string &what)
        const override;

    virtual int run() override;

    int useLexicalCast();

    int useStream();

    LexicalCast::optional lexicalCast_;
};

void MathIo::configuration(po::options_description &cmdline
                           , po::options_description &config
                           , po::positional_options_description &pd)
{
    cmdline.add_options()
        ("type", po::value<Type>()
         , "Type when using boost::lexical_cast. "
         "Must be used together with --value.")
        ("value", po::value<std::string>()
         , "Value when useing boost::lexical_cast. "
         "Must be used together with --type.")
        ;

    pd
        .add("type", 1)
        .add("value", 1);

    (void) config;
}

void MathIo::configure(const po::variables_map &vars)
{
    auto type(vars.count("type"));
    auto value(vars.count("value"));
    if (type != value) {
        throw po::error("type and value must be used together");
    }
    if (type) {
        lexicalCast_.emplace(vars["type"].as<Type>()
                             , vars["value"].as<std::string>());
    }
}

bool MathIo::help(std::ostream &out, const std::string &what) const
{
    if (what.empty()) {
        out << R"RAW(math.io: test I/O operations on various libmath data types
)RAW";
    }
    return false;
}

template <typename T>
bool parse(std::istream &is, Type type)
{
    T value;
    if (!(is >> value)) {
        std::cout << "partial parse: <" << value << ">" << std::endl;
        return false;
    }

    std::cout << "parsed " << type << "(" << value << ")" << std::endl;

    return true;
}

int MathIo::useLexicalCast()
{
    // TODO: implement me
    return EXIT_SUCCESS;
}

void state(const std::string &message, std::istream &is)
{
    std::cout << message;
    if (is.good()) { std::cout << " good"; }
    if (is.fail()) { std::cout << " fail"; }
    if (is.bad()) { std::cout << " bad"; }
    if (is.eof()) { std::cout << " eof"; }
    std::cout << std::endl;
}

int MathIo::useStream()
{
    auto &is(std::cin);

    Type type = {};
    state("> start state:", is);

    while (is >> type) {
#define CASE(TYPE)                                                      \
        case Type::TYPE:                                                \
            if (!parse<math::TYPE>(is, type)) {                         \
                std::cerr << "failed to parse " << type << std::endl;   \
                continue;                                               \
            }                                                           \
            break

        switch (type) {
            CASE(Size2);
            CASE(Size2i);
            CASE(Size2ll);
            CASE(Size2f);
            CASE(Size2r);

            CASE(Size3);
            CASE(Size3i);
            CASE(Size3ll);
            CASE(Size3f);
            CASE(Size3r);

            CASE(Extents2);
            CASE(Extents2i);
            CASE(Extents2ll);
            CASE(Extents2f);

            CASE(Extents3);
            CASE(Extents3i);
            CASE(Extents3ll);
            CASE(Extents3f);

            CASE(Viewport2);
            CASE(Viewport2i);
            CASE(Viewport2f);
        }
#undef CASE

        state("> after parse state:", is);
    }

    state("> end state:", is);

    return (is ? EXIT_SUCCESS : EXIT_FAILURE);
}

int MathIo::run()
{
    if (lexicalCast_) {
        return useLexicalCast();
    }

    return useStream();
}

} // namespace

int main(int argc, char *argv[])
{
    return MathIo()(argc, argv);
}

